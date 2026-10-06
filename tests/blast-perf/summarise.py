#!/usr/bin/env python3
"""Summarise batch-vs-per-query MEGABLAST results into REPORT.md."""
import csv
import json
import re
import sys
from pathlib import Path

NODE_COST_PER_HOUR = 7.50
EXPECTED_MAX_CONCURRENT_SINGLE = 3
MARGIN_INCONCLUSIVE_PCT = 15.0


def read_timings(outdir: Path) -> dict:
    path = outdir / "timings.json"
    if not path.exists():
        return {}
    return json.loads(path.read_text())


def read_trace(outdir: Path) -> list[dict]:
    path = outdir / "trace.txt"
    if not path.exists():
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_timing_files(outdir: Path, pattern: str) -> dict[str, dict]:
    results = {}
    for path in sorted(outdir.glob(pattern)):
        values = {}
        for line in path.read_text().splitlines():
            key, _, value = line.partition("=")
            values[key] = value
        results[path.stem] = values
    return results


def read_node_hostnames(outdir: Path) -> set[str]:
    hostnames = set()
    for path in outdir.glob("*.node.txt"):
        lines = path.read_text().splitlines()
        if lines:
            hostnames.add(lines[0])
    return hostnames


def blast_span_seconds(timings: dict[str, dict]) -> float | None:
    starts = [int(v["start"]) for v in timings.values() if v.get("start")]
    ends = [int(v["end"]) for v in timings.values() if v.get("end")]
    if not starts or not ends:
        return None
    return max(ends) - min(starts)


def sum_cpu_seconds(timings: dict[str, dict], cpus: int) -> float:
    total = 0.0
    for values in timings.values():
        seconds = values.get("seconds")
        if seconds:
            total += float(seconds) * cpus
    return total


def node_up_duration_seconds(trace_rows: list[dict]) -> float | None:
    for row in trace_rows:
        if row.get("name", "").split()[0] == "NODE_UP":
            return parse_nextflow_duration(row.get("duration"))
    return None


def max_concurrent_intervals(trace_rows: list[dict], name: str) -> int:
    intervals = []
    for row in trace_rows:
        if row.get("name", "").split()[0] != name:
            continue
        start = row.get("start")
        complete = row.get("complete")
        if (
            not start
            or not complete
            or start in ("-", "")
            or complete in ("-", "")
        ):
            continue
        intervals.append((parse_trace_time(start), parse_trace_time(complete)))

    events = []
    for start, end in intervals:
        events.append((start, 1))
        events.append((end, -1))
    events.sort()

    max_concurrent = 0
    current = 0
    for _, delta in events:
        current += delta
        max_concurrent = max(max_concurrent, current)
    return max_concurrent


def parse_trace_time(value: str) -> float:
    # trace timestamps look like "2026-10-02 07:38:33.123"; sort order is
    # all that matters here, so compare the strings' millisecond epoch
    # via fromisoformat on the normalised form.
    from datetime import datetime

    return datetime.fromisoformat(value.replace(" ", "T")).timestamp()


def diff_bls(batch_dir: Path, per_query_dir: Path) -> tuple[bool, int]:
    batch_bls = batch_dir / "batch.bls"
    if not batch_bls.exists():
        return (False, -1)

    batch_rows = sorted(batch_bls.read_text().splitlines())

    per_query_rows = []
    for path in sorted(per_query_dir.glob("queries.*.bls")):
        per_query_rows.extend(path.read_text().splitlines())
    per_query_rows.sort()

    if batch_rows == per_query_rows:
        return (True, 0)

    batch_set = set(batch_rows)
    per_query_set = set(per_query_rows)
    differing = len(batch_set ^ per_query_set)
    return (False, differing)


def parse_nextflow_duration(value: str | None) -> float | None:
    """Parse a Nextflow Duration.toString(), e.g. "1h 5m 2s" or "500ms"."""
    if not value:
        return None
    units = {"d": 86400, "h": 3600, "m": 60, "s": 1, "ms": 0.001}
    total = 0.0
    for amount, unit in re.findall(r"(\d+(?:\.\d+)?)\s*(ms|[dhms])", value):
        total += float(amount) * units[unit]
    return total or None


def format_duration(seconds: float | None) -> str:
    if seconds is None:
        return "n/a"
    minutes, secs = divmod(int(seconds), 60)
    return f"{minutes}m{secs:02d}s"


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit("usage: summarise.py <batch_outdir> <per_query_outdir>")

    batch_dir = Path(sys.argv[1])
    per_query_dir = Path(sys.argv[2])

    batch_timings = read_timings(batch_dir)
    per_query_timings = read_timings(per_query_dir)

    batch_trace = read_trace(batch_dir)
    per_query_trace = read_trace(per_query_dir)

    batch_blast_timing = read_timing_files(batch_dir, "batch.timing")
    per_query_blast_timing = read_timing_files(
        per_query_dir, "queries.*.timing"
    )

    batch_span = blast_span_seconds(batch_blast_timing)
    per_query_span = blast_span_seconds(per_query_blast_timing)

    batch_cpu_seconds = sum_cpu_seconds(batch_blast_timing, cpus=48)
    per_query_cpu_seconds = sum_cpu_seconds(per_query_blast_timing, cpus=16)

    max_concurrent = max_concurrent_intervals(per_query_trace, "BLAST_SINGLE")

    batch_hostnames = read_node_hostnames(batch_dir)
    per_query_hostnames = read_node_hostnames(per_query_dir)

    batch_provisioning = node_up_duration_seconds(batch_trace)
    per_query_provisioning = node_up_duration_seconds(per_query_trace)

    identical, differing_rows = diff_bls(batch_dir, per_query_dir)

    batch_walltime = batch_timings.get("duration")
    per_query_walltime = per_query_timings.get("duration")

    lines = []
    lines.append("# BLAST batch-vs-per-query performance report")
    lines.append("")
    lines.append(
        "**A single run of each workflow.** There is no measure of "
        "run-to-run spread, so treat any margin under "
        f"{MARGIN_INCONCLUSIVE_PCT:.0f}% as inconclusive."
    )
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append("| Metric | `batch.nf` | `per_query.nf` |")
    lines.append("|---|---|---|")
    lines.append(
        f"| Workflow walltime | {batch_walltime or 'n/a'} | "
        f"{per_query_walltime or 'n/a'} |"
    )
    lines.append(
        f"| BLAST span | {format_duration(batch_span)} | "
        f"{format_duration(per_query_span)} |"
    )
    lines.append(
        f"| Sum of blastn CPU-seconds | {batch_cpu_seconds:.0f} | "
        f"{per_query_cpu_seconds:.0f} |"
    )
    lines.append("")

    lines.append("## Per-query timing (`per_query.nf`)")
    lines.append("")
    lines.append("| Query | Seconds |")
    lines.append("|---|---|")
    for name, values in sorted(per_query_blast_timing.items()):
        lines.append(f"| {name} | {values.get('seconds', 'n/a')} |")
    lines.append("")

    concurrency_note = (
        "as expected"
        if max_concurrent == EXPECTED_MAX_CONCURRENT_SINGLE
        else f"UNEXPECTED, expected {EXPECTED_MAX_CONCURRENT_SINGLE}"
    )

    lines.append("## Checks")
    lines.append("")
    lines.append(
        f"- **Concurrency.** Max concurrent `BLAST_SINGLE` tasks: "
        f"{max_concurrent} ({concurrency_note})."
    )
    if identical:
        lines.append(
            "- **Correctness.** `batch.bls` and the concatenated "
            "per-query `.bls` files are identical."
        )
    elif differing_rows < 0:
        lines.append(
            "- **Correctness.** `batch.bls` not found; correctness "
            "check skipped."
        )
    else:
        lines.append(
            f"- **Correctness.** {differing_rows} differing rows between "
            "`batch.bls` and the concatenated per-query `.bls` files "
            "(megablast batching can shift hits at the max_target_seqs "
            "cut-off; a small difference is expected, not a failure)."
        )
    lines.append("")

    lines.append("## Node provisioning")
    lines.append("")
    lines.append(
        "Each workflow provisioned its own node (no shared warm node "
        "between runs), so they're expected to land on different "
        "hosts and pay their own staging cost."
    )
    lines.append("")
    lines.append("| | `batch.nf` | `per_query.nf` |")
    lines.append("|---|---|---|")
    lines.append(
        f"| NODE_UP duration | {format_duration(batch_provisioning)} | "
        f"{format_duration(per_query_provisioning)} |"
    )
    lines.append(
        f"| Hostname | {', '.join(sorted(batch_hostnames)) or 'n/a'} | "
        f"{', '.join(sorted(per_query_hostnames)) or 'n/a'} |"
    )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    if (
        batch_span is not None
        and per_query_span is not None
        and per_query_span > 0
    ):
        margin_pct = 100.0 * abs(batch_span - per_query_span) / per_query_span
        faster = "batch.nf" if batch_span < per_query_span else "per_query.nf"
        lines.append(
            f"**{faster}** has the shorter BLAST span on this "
            f"infrastructure, by {margin_pct:.1f}%"
            + (
                " (this margin is small enough to be inconclusive from a "
                "single run)."
                if margin_pct < MARGIN_INCONCLUSIVE_PCT
                else "."
            )
        )
    else:
        lines.append("Not enough data to compute a verdict.")
    lines.append("")

    lines.append("## Estimated cost")
    lines.append("")
    for label, duration in (
        ("batch.nf", batch_timings.get("duration")),
        ("per_query.nf", per_query_timings.get("duration")),
    ):
        seconds = parse_nextflow_duration(duration)
        cost = (
            f"${seconds / 3600 * NODE_COST_PER_HOUR:.2f}"
            if seconds is not None
            else "n/a"
        )
        lines.append(f"- {label}: {duration or 'n/a'} -> {cost}")
    lines.append(f"- Node rate: ${NODE_COST_PER_HOUR:.2f}/h")
    lines.append("")

    Path("REPORT.md").write_text("\n".join(lines) + "\n")
    print("Wrote REPORT.md")


if __name__ == "__main__":
    main()
