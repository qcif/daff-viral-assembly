#!/usr/bin/env bash
# Cloudgene post-run step: bundle the run QC report and each sample's HTML
# reports into <outdir>/reports.zip, keeping the <sample>/<name>.html layout
# so the links inside each sample report (QC report, FastQC pages) still work.
# BAM/BAI and JSON files from 11_report/ are deliberately left out.

set -euo pipefail

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

# Per-sample report dirs; skip the numbered pipeline dirs (01_*, 02_*, ...)
sample_dirs=()
for report in "${DIR}"/*/11_report/*_report.html; do
    [[ -f "$report" ]] || continue
    sample_dir="$(basename "$(dirname "$(dirname "$report")")")"
    [[ ! "$sample_dir" =~ ^0[0-9]_ ]] || continue
    if [[ " ${sample_dirs[*]:-} " != *" ${sample_dir} "* ]]; then
        sample_dirs+=("$sample_dir")
    fi
done

# Never fail this step: missing reports are flagged by collect_errors.sh.
if [[ "${#sample_dirs[@]}" -eq 0 ]]; then
    echo "WARNING: no sample reports found under ${DIR}/*/11_report/"
    echo "Nothing to zip."
    exit 0
fi

staging="$(mktemp -d)"
trap 'rm -rf "$staging"' EXIT
mkdir "${staging}/reports"

for qc in run_qc_report.html run_qc_report.txt; do
    [[ ! -f "${DIR}/02_qc_report/${qc}" ]] \
        || cp "${DIR}/02_qc_report/${qc}" "${staging}/reports/"
done

for sample_dir in "${sample_dirs[@]}"; do
    mkdir "${staging}/reports/${sample_dir}"
    for html in "${DIR}/${sample_dir}"/11_report/*.html; do
        [[ -f "$html" ]] || continue
        cp "$html" "${staging}/reports/${sample_dir}/"
    done
done

rm -f "${DIR}/reports.zip"
(cd "$staging" && zip -qrX "$(realpath "$DIR")/reports.zip" reports)

echo "All workflow reports have been zipped"
echo 'Download all reports from reports.zip under the "Results" tab.'
