#!/bin/bash
# Fetches ~30,000 viral nucleotide sequences, 9-11kb, for the scale test
# (scale_batch.nf / scale_split.nf): one big batched blastn call vs. many
# smaller parallel ones, on realistic query volume instead of the 10
# queries used in batch-vs-series.md and cores-allocated.md.
#
# Output is NOT committed (it's ~250-350MB): everything lands under
# scale_queries/, which .gitignore excludes. Re-run this script to
# regenerate it; the exact records returned can drift slightly between
# runs as NCBI's database changes, but that doesn't matter for a
# throughput/scaling test the way it would for a correctness test.
#
# No NCBI API key is configured for this project, so requests are capped
# at the unauthenticated rate limit (3/sec) via SLEEP_BETWEEN below.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_DIR="${SCRIPT_DIR}/scale_queries"
OUT_FASTA="${OUT_DIR}/queries_30k.fasta"
OUT_COUNT="${OUT_DIR}/queries_30k.n"
EUTILS_BASE_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
TARGET_COUNT=30000
BATCH_SIZE=500
MIN_LEN=9000
MAX_LEN=11000
SLEEP_BETWEEN=0.4
# Genomic (not mRNA) viral sequences in the target length range. This
# pulls from all of nuccore, not just RefSeq (RefSeq alone only has ~900
# records in this range, nowhere near enough for 30,000).
SEARCH_TERM='Viruses[Organism] AND 9000:11000[SLEN] AND biomol_genomic[PROP]'

mkdir -p "${OUT_DIR}"

echo "Searching nuccore for: ${SEARCH_TERM}" >&2
search_xml="$(curl -s --max-time 30 \
    --data-urlencode "db=nuccore" \
    --data-urlencode "term=${SEARCH_TERM}" \
    --data-urlencode "usehistory=y" \
    --data-urlencode "retmax=0" \
    "${EUTILS_BASE_URL}/esearch.fcgi")"

count=$(grep -oP '(?<=<Count>)\d+' <<< "${search_xml}" | head -1)
webenv=$(grep -oP '(?<=<WebEnv>)[^<]+' <<< "${search_xml}")
query_key=$(grep -oP '(?<=<QueryKey>)[^<]+' <<< "${search_xml}")

if [[ -z "${count}" || -z "${webenv}" || -z "${query_key}" ]]; then
    echo "ERROR: esearch didn't return Count/WebEnv/QueryKey:" >&2
    echo "${search_xml}" >&2
    exit 1
fi

echo "Found ${count} matching records; fetching up to ${TARGET_COUNT}" >&2

fetch_count=$((count < TARGET_COUNT ? count : TARGET_COUNT))
: > "${OUT_FASTA}"

retstart=0
while (( retstart < fetch_count )); do
    retmax=$BATCH_SIZE
    if (( retstart + retmax > fetch_count )); then
        retmax=$((fetch_count - retstart))
    fi
    echo "Fetching records ${retstart}-$((retstart + retmax))..." >&2
    curl -s --max-time 60 \
        --data-urlencode "db=nuccore" \
        --data-urlencode "WebEnv=${webenv}" \
        --data-urlencode "query_key=${query_key}" \
        --data-urlencode "retstart=${retstart}" \
        --data-urlencode "retmax=${retmax}" \
        --data-urlencode "rettype=fasta" \
        --data-urlencode "retmode=text" \
        "${EUTILS_BASE_URL}/efetch.fcgi" >> "${OUT_FASTA}"
    retstart=$((retstart + retmax))
    sleep "${SLEEP_BETWEEN}"
done

# Drop any record outside the length window (SLEN is indexed on the full
# record; a handful of edge cases can still slip through) and any
# truncated/empty record from a dropped connection.
echo "Filtering to ${MIN_LEN}-${MAX_LEN}bp and de-duplicating..." >&2
../../venv/bin/python3 - "${OUT_FASTA}" "${MIN_LEN}" "${MAX_LEN}" <<'PYEOF'
import sys

path, min_len, max_len = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
records = []
header, seq_chunks = None, []
seen_ids = set()

with open(path) as fh:
    for line in fh:
        if line.startswith(">"):
            if header is not None:
                records.append((header, "".join(seq_chunks)))
            header, seq_chunks = line.rstrip("\n"), []
        else:
            seq_chunks.append(line.strip())
    if header is not None:
        records.append((header, "".join(seq_chunks)))

kept = []
for header, seq in records:
    acc = header.split()[0][1:]
    if acc in seen_ids:
        continue
    if not (min_len <= len(seq) <= max_len):
        continue
    seen_ids.add(acc)
    kept.append((header, seq))

with open(path, "w") as fh:
    for header, seq in kept:
        fh.write(header + "\n")
        for i in range(0, len(seq), 70):
            fh.write(seq[i:i + 70] + "\n")

print(f"Kept {len(kept)}/{len(records)} records", file=sys.stderr)
with open(sys.argv[1].rsplit('/', 1)[0] + "/queries_30k.n", "w") as fh:
    fh.write(str(len(kept)))
PYEOF

n=$(cat "${OUT_COUNT}")
echo "Wrote ${OUT_FASTA} (${n} sequences) and ${OUT_COUNT}" >&2
