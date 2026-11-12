#!/usr/bin/env bash
# Cloudgene post-run step: summarise task and sample status from the
# Nextflow execution trace and exit non-zero if any task failed or any
# sample has no report.

set -euo pipefail

EXIT_NO_TRACE=3
EXIT_PROBLEMS=1

usage() {
    echo "Usage: $0 --dir <outdir> [--samplesheet <samplesheet.csv>]"
    echo "Reads <outdir>/01_pipeline_logs/execution_trace_*.txt"
}

DIR=""
SAMPLESHEET=""
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --dir)
            [[ -n "${2:-}" ]] || { echo "Missing value for --dir"; usage; exit 1; }
            DIR="${2%/}"
            shift 2
            ;;
        --samplesheet)
            [[ -n "${2:-}" ]] || { echo "Missing value for --samplesheet"; usage; exit 1; }
            SAMPLESHEET="$2"
            shift 2
            ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
done
[[ -n "$DIR" ]] || { usage; exit 1; }

logs_dir="${DIR}/01_pipeline_logs"

# Timestamped names (YYYY-MM-DD_HH-MM-SS) sort lexically; take the newest.
trace=""
if [[ -d "$logs_dir" ]]; then
    trace="$(find "$logs_dir" -maxdepth 1 -type f \
        -name 'execution_trace_*.txt' | sort | tail -n 1)"
fi

if [[ -z "$trace" || ! -s "$trace" ]]; then
    echo "ERROR: no execution trace found in ${logs_dir}."
    echo "The run did not finish writing its results; see the workflow log."
    exit "$EXIT_NO_TRACE"
fi

# Task status counts, failed tasks, then per-sample report status.
# Exits 1 if any task failed or any sample lacks a report.
rc=0
python3 - "$DIR" "$trace" "$SAMPLESHEET" <<'PY' || rc=$?
import csv
import os
import re
import sys
from collections import Counter

BAD_STATUSES = ('FAILED', 'ABORTED')
OK_STATUSES = ('COMPLETED', 'CACHED')
WORK_BASE = 'az://workdata/wf4'
NATIVE_ID_RE = re.compile(r'nf-([0-9a-f]{32})\s*$')

outdir, trace, samplesheet = sys.argv[1:4]
problems = False

with open(trace, newline='') as handle:
    tasks = list(csv.DictReader(handle, delimiter='\t'))

print(f'Task status summary (from {os.path.basename(trace)})')
counts = Counter(t.get('status', '') for t in tasks)
for status, count in sorted(counts.items()):
    print(f'  {status}: {count}')

# A failed attempt that was retried successfully keeps its FAILED row.
succeeded = {
    t.get('name') for t in tasks
    if t.get('status') in OK_STATUSES
}

for task in tasks:
    if task.get('status') not in BAD_STATUSES:
        continue
    retried = task.get('name') in succeeded
    if task['status'] == 'FAILED' and not retried:
        problems = True
    print(f"\n{task.get('name')} [{task['status']}]: "
          f"exit {task.get('exit')}"
          + (' (retried successfully)' if retried else ''))
    match = NATIVE_ID_RE.search(task.get('native_id') or '')
    if match:
        digest = match.group(1)
        print(f'  work dir: {WORK_BASE}/{digest[:2]}/{digest[2:]}')

# Sample list: samplesheet if readable, else dirs holding an 11_report/.
samples = []
if samplesheet:
    try:
        with open(samplesheet, newline='') as handle:
            for row in csv.DictReader(handle):
                sample = (row.get('sample') or '').strip()
                if sample and sample not in samples:
                    samples.append(sample)
    except OSError:
        print(f'\nWARNING: cannot read samplesheet {samplesheet}; '
              'listing samples from the output directory')
if not samples:
    for name in sorted(os.listdir(outdir)):
        if os.path.isdir(os.path.join(outdir, name, '11_report')):
            samples.append(name)

flags = {}
qc_path = os.path.join(outdir, '02_qc_report', 'run_qc_report.txt')
if os.path.isfile(qc_path):
    with open(qc_path, newline='') as handle:
        for row in csv.DictReader(handle, delimiter='\t'):
            if row.get('Sample') and row.get('QC_FLAG'):
                flags[row['Sample']] = row['QC_FLAG']

print('\nSample status summary')
for sample in samples:
    report = os.path.join(
        outdir, sample, '11_report', f'{sample}_report.html')
    if os.path.isfile(report) and os.path.getsize(report) > 0:
        line = f'  {sample}: report OK'
    else:
        line = f'  {sample}: NO REPORT'
        problems = True
    if sample in flags:
        line += f'  QC: {flags[sample]}'
    print(line)

sys.exit(1 if problems else 0)
PY

echo
if [[ "$rc" -eq 0 ]]; then
    echo "No task errors"
    exit 0
elif [[ "$rc" -eq "$EXIT_PROBLEMS" ]]; then
    echo "Check the workflow step log for details of the failures above."
    exit "$EXIT_PROBLEMS"
fi
echo "ERROR: failed to parse ${trace}; check the workflow step log."
exit "$rc"
