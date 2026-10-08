#!/usr/bin/env python3
"""Summarise the thread-count sweep (cores.nf) into REPORT.md."""
import sys
from pathlib import Path

from summarise import (
    format_duration,
    node_up_duration_seconds,
    parse_nextflow_duration,
    read_node_hostnames,
    read_timing_files,
    read_timings,
    read_trace,
)

NODE_COST_PER_HOUR = 7.50
NODE_CPUS = 48
CLOCK_TICKS_PER_SEC = 100.0


def read_samples(path: Path) -> list[tuple[int, float, int]]:
    """Read a `<pass>.sample` file: (unix_time, cpu_ticks, disk_sectors)."""
    samples = []
    if not path.exists():
        return samples
    for line in path.read_text().splitlines():
        parts = line.split()
        if len(parts) != 3:
            continue
        samples.append((int(parts[0]), float(parts[1]), int(parts[2])))
    return samples


def peak_cpu_seconds_per_sec(samples: list[tuple[int, float, int]]) -> float:
    """Max instantaneous CPU-seconds/sec between consecutive samples."""
    peak = 0.0
    for (t0, cpu0, _), (t1, cpu1, _) in zip(samples, samples[1:]):
        dt = t1 - t0
        if dt <= 0:
            continue
        rate = (cpu1 - cpu0) / CLOCK_TICKS_PER_SEC / dt
        peak = max(peak, rate)
    return peak


def disk_sectors_read(samples: list[tuple[int, float, int]]) -> int:
    if len(samples) < 2:
        return 0
    return samples[-1][2] - samples[0][2]


def total_cpu_seconds(samples: list[tuple[int, float, int]]) -> float:
    if len(samples) < 2:
        return 0.0
    return (samples[-1][1] - samples[0][1]) / CLOCK_TICKS_PER_SEC


def read_cgroup_limits(outdir: Path) -> dict[str, str]:
    """cpu.max / memory.max are the same across passes (one node); take the
    first pass's node.txt. Its last two lines are always `cpu.max`'s value
    (or `cpu.max=n/a`) followed by `memory.max`'s value (or
    `memory.max=n/a`) — see the `BLAST_PASS` script in modules.nf."""
    limits = {"cpu.max": "n/a", "memory.max": "n/a"}
    paths = sorted(outdir.glob("*.node.txt"))
    for path in paths:
        if path.name == "node_up.node.txt":
            continue
        lines = [ln for ln in path.read_text().splitlines() if ln]
        if len(lines) < 2:
            continue
        cpu_max, mem_max = lines[-2], lines[-1]
        limits["cpu.max"] = (
            "n/a" if cpu_max == "cpu.max=n/a" else cpu_max
        )
        limits["memory.max"] = (
            "n/a" if mem_max == "memory.max=n/a" else mem_max
        )
        break
    return limits


def diff_all_bls(outdir: Path) -> dict[str, int]:
    """Diff every pass's .bls against cold_48.bls. Returns {name: differing
    row count}, 0 meaning identical."""
    cold_path = outdir / "cold_48.bls"
    if not cold_path.exists():
        return {}
    cold_rows = set(cold_path.read_text().splitlines())
    results = {}
    for path in sorted(outdir.glob("warm_*.bls")):
        rows = set(path.read_text().splitlines())
        results[path.stem] = len(cold_rows ^ rows)
    return results


def main() -> None:
    if len(sys.argv) != 2:
        sys.exit("usage: summarise_cores.py <cores_outdir>")

    outdir = Path(sys.argv[1])

    timings = read_timings(outdir)
    trace = read_trace(outdir)
    blast_timing = read_timing_files(outdir, "*_*.timing")

    def sort_key(item: tuple[str, dict]) -> tuple[int, int]:
        name, values = item
        kind_rank = 0 if name.startswith("cold") else 1
        threads = values.get("threads", "0")
        threads_n = int(threads) if threads.isdigit() else 0
        return (kind_rank, threads_n)

    passes = []
    for name, values in sorted(blast_timing.items(), key=sort_key):
        threads = values.get("threads", "?")
        kind = "cold" if name.startswith("cold") else "warm"
        samples = read_samples(outdir / f"{name}.sample")
        passes.append({
            "name": name,
            "kind": kind,
            "threads": threads,
            "seconds": values.get("seconds"),
            "cpu_seconds": total_cpu_seconds(samples),
            "peak_cpu_rate": peak_cpu_seconds_per_sec(samples),
            "disk_sectors": disk_sectors_read(samples),
        })

    hostnames = read_node_hostnames(outdir)
    provisioning = node_up_duration_seconds(trace)
    cgroup_limits = read_cgroup_limits(outdir)
    bls_diff = diff_all_bls(outdir)
    walltime = timings.get("duration")

    lines = []
    lines.append("# BLAST thread-count sweep report")
    lines.append("")
    lines.append(
        "**One cold pass and one run of each warm thread count.** There "
        "is no measure of run-to-run spread, so treat any difference "
        "between adjacent thread counts as inconclusive."
    )
    lines.append("")
    lines.append("## Sweep")
    lines.append("")
    lines.append(
        "| Pass | Threads | blastn time | CPU-seconds used | "
        "Peak CPU-sec/sec | Disk sectors read |"
    )
    lines.append("|---|---|---|---|---|---|")
    for p in passes:
        seconds = p["seconds"]
        duration = format_duration(float(seconds)) if seconds else "n/a"
        lines.append(
            f"| {p['kind']} | {p['threads']} | {duration} | "
            f"{p['cpu_seconds']:.0f} | {p['peak_cpu_rate']:.1f} | "
            f"{p['disk_sectors']} |"
        )
    lines.append("")

    cold_passes = [p for p in passes if p["kind"] == "cold"]
    warm_passes = [p for p in passes if p["kind"] == "warm"]

    lines.append("## Checks")
    lines.append("")
    if bls_diff:
        all_identical = all(v == 0 for v in bls_diff.values())
        if all_identical:
            lines.append(
                "- **Correctness.** Every warm pass's `.bls` is identical "
                "to the cold pass's. Thread count doesn't change hits."
            )
        else:
            for name, count in bls_diff.items():
                note = "identical" if count == 0 else f"{count} differing rows"
                lines.append(f"- **Correctness ({name}).** {note}.")
    else:
        lines.append("- **Correctness.** `cold_48.bls` not found; skipped.")

    if cold_passes and warm_passes:
        cold_disk = cold_passes[0]["disk_sectors"]
        suspect = [
            p for p in warm_passes
            if cold_disk > 0 and p["disk_sectors"] > 0.5 * cold_disk
        ]
        if suspect:
            names = ", ".join(p["name"] for p in suspect)
            lines.append(
                f"- **Warm really was warm?** {names} read close to as "
                "many disk sectors as the cold pass — cache eviction may "
                "have happened mid-sweep; treat those points as suspect."
            )
        else:
            lines.append(
                "- **Warm really was warm.** Every warm pass's disk reads "
                "were well below the cold pass's."
            )
    lines.append("")

    lines.append("## Node")
    lines.append("")
    lines.append(
        f"- NODE_UP duration: {format_duration(provisioning)}, "
        f"hostname: {', '.join(sorted(hostnames)) or 'n/a'}"
    )
    lines.append(f"- `cpu.max`: `{cgroup_limits['cpu.max']}`")
    lines.append(f"- `memory.max`: `{cgroup_limits['memory.max']}`")
    lines.append(f"- Workflow walltime: {walltime or 'n/a'}")
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    if len(warm_passes) >= 2:
        baseline = warm_passes[0]
        baseline_seconds = (
            float(baseline["seconds"]) if baseline["seconds"] else None
        )
        if baseline_seconds:
            lines.append(
                f"Speedup relative to {baseline['threads']} threads "
                "(warm):"
            )
            lines.append("")
            lines.append(
                "| Threads | blastn time | Speedup | Tasks/node (warm) |"
            )
            lines.append("|---|---|---|---|")
            for p in warm_passes:
                seconds = float(p["seconds"]) if p["seconds"] else None
                speedup = (
                    f"{baseline_seconds / seconds:.2f}x"
                    if seconds else "n/a"
                )
                threads_int = (
                    int(p["threads"]) if str(p["threads"]).isdigit() else None
                )
                tasks_per_node = (
                    f"{NODE_CPUS // threads_int}" if threads_int else "n/a"
                )
                lines.append(
                    f"| {p['threads']} | "
                    f"{format_duration(seconds) if seconds else 'n/a'} | "
                    f"{speedup} | {tasks_per_node} |"
                )
            lines.append("")
            lines.append(
                "\"Tasks/node\" assumes CPU count alone determines how "
                "many tasks of that thread count fit on a 48-core node "
                "concurrently (ignoring memory); multiply blastn time's "
                "inverse by that to compare total throughput per node."
            )
        else:
            lines.append("Not enough data to compute a verdict.")
    else:
        lines.append("Not enough data to compute a verdict.")
    lines.append("")

    lines.append("## Estimated cost")
    lines.append("")
    seconds = parse_nextflow_duration(walltime)
    cost = (
        f"${seconds / 3600 * NODE_COST_PER_HOUR:.2f}"
        if seconds is not None
        else "n/a"
    )
    lines.append(f"- Workflow: {walltime or 'n/a'} -> {cost}")
    lines.append(f"- Node rate: ${NODE_COST_PER_HOUR:.2f}/h")
    lines.append("")

    Path("REPORT.md").write_text("\n".join(lines) + "\n")
    print("Wrote REPORT.md")


if __name__ == "__main__":
    main()
