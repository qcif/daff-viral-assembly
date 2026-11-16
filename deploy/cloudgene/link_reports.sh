#!/usr/bin/env bash
# Cloudgene post-run step: gather each sample's HTML report into a single
# <outdir>/reports/ folder so users can find them without opening every
# sample directory.
#
# Sample reports link to sibling pages by relative path (./run_qc_report.html,
# ./<sample>_*_fastqc.html), so those pages are linked into reports/ too. All
# of them are sample-prefixed apart from the run QC report, which is shared,
# so a flat folder keeps every link working.
#
# Files are hard-linked (same filesystem, no extra space) so they are served
# as ordinary files; falls back to copying if hard links are not possible.

set -euo pipefail

REPORTS_DIRNAME="reports"
RUN_QC_REPORT="02_qc_report/run_qc_report.html"

usage() {
    echo "Usage: $0 --dir <outdir>"
}

DIR=""
while [[ "$#" -gt 0 ]]; do
    case "$1" in
        --dir)
            [[ -n "${2:-}" ]] || { echo "Missing value for --dir"; usage; exit 1; }
            DIR="${2%/}"
            shift 2
            ;;
        -h|--help) usage; exit 0 ;;
        *) echo "Unknown parameter passed: $1"; usage; exit 1 ;;
    esac
done
[[ -n "$DIR" ]] || { usage; exit 1; }

reports_dir="${DIR}/${REPORTS_DIRNAME}"

link_file() {
    ln -f "$1" "${reports_dir}/" 2>/dev/null || cp -f "$1" "${reports_dir}/"
}

# Per-sample report dirs; skip the numbered pipeline dirs (01_*, 02_*, ...)
sample_dirs=()
for report in "${DIR}"/*/11_report/*_report.html; do
    [[ -f "$report" ]] || continue
    sample_dir="$(dirname "$(dirname "$report")")"
    [[ ! "$(basename "$sample_dir")" =~ ^0[0-9]_ ]] || continue
    if [[ " ${sample_dirs[*]:-} " != *" ${sample_dir} "* ]]; then
        sample_dirs+=("$sample_dir")
    fi
done

# Never fail this step: missing reports are flagged by collect_errors.sh.
if [[ "${#sample_dirs[@]}" -eq 0 ]]; then
    echo "WARNING: no sample reports found under ${DIR}/*/11_report/"
    echo "Nothing to link."
    exit 0
fi

mkdir -p "$reports_dir"

[[ ! -f "${DIR}/${RUN_QC_REPORT}" ]] || link_file "${DIR}/${RUN_QC_REPORT}"

for sample_dir in "${sample_dirs[@]}"; do
    for html in "${sample_dir}"/11_report/*.html; do
        [[ -f "$html" ]] || continue
        # Per-sample copy of the run QC report; already linked above.
        if [[ "$(basename "$html")" == "$(basename "$RUN_QC_REPORT")" \
              && -f "${reports_dir}/$(basename "$RUN_QC_REPORT")" ]]; then
            continue
        fi
        link_file "$html"
    done
done

echo "Linked ${#sample_dirs[@]} sample report(s) into ${REPORTS_DIRNAME}/"
echo "Open ${REPORTS_DIRNAME}/<sample>_report.html under the \"Results\" tab."
