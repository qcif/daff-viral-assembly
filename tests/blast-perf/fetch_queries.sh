#!/bin/bash
# Fetches the RefSeq genomes used as BLAST queries for the batch-vs-per-query
# performance test, rewrites headers to short sortable IDs, and writes
# queries.tsv. Regenerate queries.fasta with this script rather than editing
# it by hand, so both test workflows keep using byte-identical input.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT_FASTA="${SCRIPT_DIR}/queries.fasta"
OUT_TSV="${SCRIPT_DIR}/queries.tsv"
MIN_LEN=9000
MAX_LEN=11500
EUTILS_BASE_URL="https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# accession:name:family
QUERIES=(
    "NC_001616:Potato virus Y:Potyviridae"
    "NC_001445:Plum pox virus:Potyviridae"
    "NC_002509:Turnip mosaic virus:Potyviridae"
    "NC_003224:Zucchini yellow mosaic virus:Potyviridae"
    "NC_002634:Soybean mosaic virus:Potyviridae"
    "NC_001886:Wheat streak mosaic virus:Potyviridae"
    "NC_001477:Dengue virus 1:Flaviviridae"
    "NC_012532:Zika virus:Flaviviridae"
    "NC_009942:West Nile virus:Flaviviridae"
    "NC_002031:Yellow fever virus:Flaviviridae"
)

fetch_fasta() {
    local accession="$1"
    if command -v efetch &>/dev/null; then
        efetch -db nucleotide -id "${accession}" -format fasta
    else
        curl -s "${EUTILS_BASE_URL}/efetch.fcgi?db=nucleotide&id=${accession}&rettype=fasta&retmode=text"
    fi
}

: > "${OUT_FASTA}"
: > "${OUT_TSV}"
echo -e "accession\tname\tfamily\tlength" >> "${OUT_TSV}"

i=0
for entry in "${QUERIES[@]}"; do
    i=$((i + 1))
    accession="$(cut -d: -f1 <<< "${entry}")"
    name="$(cut -d: -f2 <<< "${entry}")"
    family="$(cut -d: -f3 <<< "${entry}")"
    qid="$(printf 'q%02d_%s' "${i}" "${accession}")"

    echo "Fetching ${accession} (${name})..." >&2
    raw_fasta="$(fetch_fasta "${accession}")"
    if [[ -z "${raw_fasta}" ]]; then
        echo "ERROR: empty response fetching ${accession}" >&2
        exit 1
    fi

    seq="$(grep -v '^>' <<< "${raw_fasta}" | tr -d '[:space:]')"
    length=${#seq}

    if (( length < MIN_LEN || length > MAX_LEN )); then
        echo "ERROR: ${accession} length ${length}bp is outside" \
            "${MIN_LEN}-${MAX_LEN}bp; replace this accession" >&2
        exit 1
    fi

    echo ">${qid}" >> "${OUT_FASTA}"
    fold -w 70 <<< "${seq}" >> "${OUT_FASTA}"
    echo -e "${accession}\t${name}\t${family}\t${length}" >> "${OUT_TSV}"
done

echo "Wrote ${OUT_FASTA} and ${OUT_TSV}" >&2
