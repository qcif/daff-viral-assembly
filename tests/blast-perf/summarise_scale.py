#!/usr/bin/env python3
"""Summarise scale_batch-vs-scale_split results (~30k queries) into
REPORT.md."""
import sys
from pathlib import Path

from summarise import (
    blast_span_seconds,
    format_duration,
    max_concurrent_intervals,
    node_up_duration_seconds,
    parse_nextflow_duration,
    read_node_hostnames,
    read_timing_files,
    read_timings,
    read_trace,
    sum_cpu_seconds,
)

NODE_COST_PER_HOUR = 7.50
EXPECTED_MAX_CONCURRENT_CHUNKS = 3
MARGIN_INCONCLUSIVE_PCT = 15.0


def diff_bls(batch_dir: Path, split_dir: Path) -> tuple[bool, int]:
    batch_bls = batch_dir / "scale_batch.bls"
    if not batch_bls.exists():
        return (False, -1)

    batch_rows = set(batch_bls.read_text().splitlines())

    split_rows = set()
    for path in sorted(split_dir.glob("*.bls")):
        split_rows.update(path.read_text().splitlines())

    if batch_rows == split_rows:
        return (True, 0)
    return (False, len(batch_rows ^ split_rows))


def main() -> None:
    if len(sys.argv) != 3:
        sys.exit(
            "usage: summarise_scale.py <scale_batch_outdir> "
            "<scale_split_outdir>"
        )

    batch_dir = Path(sys.argv[1])
    split_dir = Path(sys.argv[2])

    batch_timings = read_timings(batch_dir)
    split_timings = read_timings(split_dir)

    batch_trace = read_trace(batch_dir)
    split_trace = read_trace(split_dir)

    batch_blast_timing = read_timing_files(batch_dir, "scale_batch.timing")
    split_blast_timing = read_timing_files(split_dir, "*.timing")

    batch_span = blast_span_seconds(batch_blast_timing)
    split_span = blast_span_seconds(split_blast_timing)

    batch_cpu_seconds = sum_cpu_seconds(batch_blast_timing, cpus=48)
    split_cpu_seconds = sum_cpu_seconds(split_blast_timing, cpus=16)

    max_concurrent = max_concurrent_intervals(split_trace, "SCALE_SPLIT_CHUNK")

    batch_hostnames = read_node_hostnames(batch_dir)
    split_hostnames = read_node_hostnames(split_dir)

    batch_provisioning = node_up_duration_seconds(batch_trace)
    split_provisioning = node_up_duration_seconds(split_trace)

    identical, differing_rows = diff_bls(batch_dir, split_dir)

    batch_walltime = batch_timings.get("duration")
    split_walltime = split_timings.get("duration")

    n_chunks = len(split_blast_timing)

    lines = []
    lines.append("# BLAST scale test report (~30k queries)")
    lines.append("")
    lines.append(
        "**A single run of each workflow.** There is no measure of "
        "run-to-run spread, so treat any margin under "
        f"{MARGIN_INCONCLUSIVE_PCT:.0f}% as inconclusive. Each workflow "
        "ran on its own fresh node, as in batch-vs-series.md."
    )
    lines.append("")
    lines.append("## Summary")
    lines.append("")
    lines.append(
        "| Metric | `scale_batch.nf` (48C x1) | `scale_split.nf` (16C x3) |"
    )
    lines.append("|---|---|---|")
    lines.append(
        f"| Workflow walltime | {batch_walltime or 'n/a'} | "
        f"{split_walltime or 'n/a'} |"
    )
    lines.append(
        f"| BLAST span | {format_duration(batch_span)} | "
        f"{format_duration(split_span)} |"
    )
    lines.append(
        f"| Sum of blastn CPU-seconds | {batch_cpu_seconds:.0f} | "
        f"{split_cpu_seconds:.0f} |"
    )
    lines.append(f"| Chunks | 1 | {n_chunks} |")
    lines.append("")

    lines.append("## Per-chunk timing (`scale_split.nf`)")
    lines.append("")
    lines.append("| Chunk | Seconds |")
    lines.append("|---|---|")
    for name, values in sorted(split_blast_timing.items()):
        lines.append(f"| {name} | {values.get('seconds', 'n/a')} |")
    lines.append("")

    concurrency_note = (
        "as expected"
        if max_concurrent == EXPECTED_MAX_CONCURRENT_CHUNKS
        else f"UNEXPECTED, expected {EXPECTED_MAX_CONCURRENT_CHUNKS}"
    )

    lines.append("## Checks")
    lines.append("")
    lines.append(
        f"- **Concurrency.** Max concurrent `SCALE_SPLIT_CHUNK` tasks: "
        f"{max_concurrent} ({concurrency_note})."
    )
    if identical:
        lines.append(
            "- **Correctness.** `scale_batch.bls` and the concatenated "
            "split-chunk `.bls` files are identical."
        )
    elif differing_rows < 0:
        lines.append(
            "- **Correctness.** `scale_batch.bls` not found; correctness "
            "check skipped."
        )
    else:
        lines.append(
            f"- **Correctness.** {differing_rows} differing rows between "
            "`scale_batch.bls` and the concatenated split-chunk `.bls` "
            "files (megablast batching can shift hits at the "
            "max_target_seqs cut-off; a small difference is expected, "
            "not a failure)."
        )
    lines.append("")

    lines.append("## Node provisioning")
    lines.append("")
    lines.append(
        "Each workflow provisioned its own node (no shared warm node "
        "between runs), so they're expected to land on different hosts "
        "and pay their own staging cost."
    )
    lines.append("")
    lines.append("| | `scale_batch.nf` | `scale_split.nf` |")
    lines.append("|---|---|---|")
    lines.append(
        f"| NODE_UP duration | {format_duration(batch_provisioning)} | "
        f"{format_duration(split_provisioning)} |"
    )
    lines.append(
        f"| Hostname | {', '.join(sorted(batch_hostnames)) or 'n/a'} | "
        f"{', '.join(sorted(split_hostnames)) or 'n/a'} |"
    )
    lines.append("")

    lines.append("## Verdict")
    lines.append("")
    if batch_span is not None and split_span is not None and split_span > 0:
        margin_pct = 100.0 * abs(batch_span - split_span) / split_span
        faster = (
            "scale_batch.nf" if batch_span < split_span else "scale_split.nf"
        )
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
        ("scale_batch.nf", batch_timings.get("duration")),
        ("scale_split.nf", split_timings.get("duration")),
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
