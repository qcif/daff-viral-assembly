#!/usr/bin/env bash
# Build a small "dev" Kaiju database (~500 MB target) for wf4.
# Idempotent: rerun to resume after failures. Delete $OUT_ROOT/kaiju_db to
# force a full rebuild.
#
# Env overrides:
#   OUT_ROOT       Output root (default: <repo>/refdata_dev). Must be writable.
#   KAIJU_IMAGE    Docker image (default set in lib.sh)

set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib.sh
source "$HERE/lib.sh"

resolve_out_root "$HERE"

SEED_FILE="$HERE/seed_taxa.txt"
DB_DIR="$OUT_ROOT/kaiju_db"
WORK_DIR="$OUT_ROOT/.build/kaiju"
FASTA_DIR="$WORK_DIR/proteins"
TAXO_DIR="$WORK_DIR/taxonomy"

VIRAL_PROTEIN_URL="https://ftp.ncbi.nlm.nih.gov/refseq/release/viral/viral.1.protein.faa.gz"
TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

FINAL_FMI="$DB_DIR/kaiju_dev_db.fmi"

require_cmd docker curl awk gzip tar
ensure_image "$KAIJU_IMAGE"

log "OUT_ROOT = $OUT_ROOT"

if [[ -f "$FINAL_FMI" && -f "$DB_DIR/nodes.dmp" && -f "$DB_DIR/names.dmp" ]]; then
    log "Kaiju dev DB already present at $DB_DIR — nothing to do."
    du -sh "$DB_DIR" >&2
    exit 0
fi

mkdir -p "$DB_DIR" "$FASTA_DIR" "$TAXO_DIR"

# ---- taxonomy -----------------------------------------------------------
if [[ ! -s "$DB_DIR/nodes.dmp" || ! -s "$DB_DIR/names.dmp" ]]; then
    log "Downloading NCBI taxdump"
    if [[ ! -s "$TAXO_DIR/taxdump.tar.gz" ]]; then
        curl -fsSL "$TAXDUMP_URL" -o "$TAXO_DIR/taxdump.tar.gz.tmp"
        mv "$TAXO_DIR/taxdump.tar.gz.tmp" "$TAXO_DIR/taxdump.tar.gz"
    fi
    tar -xzf "$TAXO_DIR/taxdump.tar.gz" -C "$TAXO_DIR" nodes.dmp names.dmp merged.dmp
    cp "$TAXO_DIR/nodes.dmp" "$TAXO_DIR/names.dmp" "$DB_DIR/"
fi

# ---- viral proteins (full viral RefSeq protein) --------------------------
VIRAL_FAA_GZ="$FASTA_DIR/viral.protein.faa.gz"
if [[ ! -s "$VIRAL_FAA_GZ" ]]; then
    log "Downloading viral RefSeq protein"
    curl -fsSL "$VIRAL_PROTEIN_URL" -o "$VIRAL_FAA_GZ.tmp"
    mv "$VIRAL_FAA_GZ.tmp" "$VIRAL_FAA_GZ"
fi

# ---- fabaceae + context proteins (efetch linked proteins per accession) --
for section in fabaceae context; do
    while IFS= read -r acc; do
        [[ -z "$acc" ]] && continue
        efetch_linked_proteins "$acc" "$FASTA_DIR/${acc}.faa"
    done < <(parse_section "$section" "$SEED_FILE")
done

# ---- Build the accession -> taxid map ------------------------------------
# Kaiju's FMI format needs each protein header rewritten to ">N_TAXID" where N
# is a unique integer per sequence. Resolve taxid per unique accession via
# NCBI esummary (cached to $ACC2TAXID).
ACC2TAXID="$WORK_DIR/acc2taxid.tsv"
touch "$ACC2TAXID"

resolve_taxid() {
    local acc="$1"
    local hit
    hit=$(awk -v a="$acc" '$1 == a {print $2; exit}' "$ACC2TAXID")
    if [[ -n "$hit" ]]; then
        echo "$hit"
        return 0
    fi
    hit=$(curl -fsSL \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=${acc}" \
        | grep -oE '<Item Name="TaxId"[^>]*>[0-9]+</Item>' \
        | grep -oE '[0-9]+' | head -n1)
    [[ -z "$hit" ]] && hit=1
    printf '%s\t%s\n' "$acc" "$hit" >> "$ACC2TAXID"
    echo "$hit"
}

MERGED_FAA="$WORK_DIR/kaiju_input.faa"
: > "$MERGED_FAA"
seq_counter=0

log "Rewriting viral RefSeq headers"
zcat "$VIRAL_FAA_GZ" | awk -v out="$MERGED_FAA" '
    BEGIN { n = 0 }
    /^>/ {
        split(substr($0, 2), a, " "); acc = a[1]
        n++
        printf(">%d__ACC__%s\n", n, acc) >> out
        next
    }
    { print >> out }
'
seq_counter=$(grep -c '^>' "$MERGED_FAA" || true)
log "  viral protein sequences: $seq_counter"

for faa in "$FASTA_DIR"/*.faa; do
    [[ -s "$faa" ]] || continue
    while IFS= read -r line; do
        if [[ $line == ">"* ]]; then
            seq_counter=$((seq_counter + 1))
            hdr_acc=$(echo "$line" | awk '{print substr($1,2)}')
            printf '>%d__ACC__%s\n' "$seq_counter" "$hdr_acc" >> "$MERGED_FAA"
        else
            printf '%s\n' "$line" >> "$MERGED_FAA"
        fi
    done < "$faa"
done

log "Resolving taxids for unique accessions (hits NCBI E-utilities, cached)"
mapfile -t uniq_accs < <(grep -oE '__ACC__[^ ]+' "$MERGED_FAA" | sed 's/__ACC__//' | sort -u)
log "  ${#uniq_accs[@]} unique accessions"
i=0
for acc in "${uniq_accs[@]}"; do
    resolve_taxid "$acc" >/dev/null
    i=$((i + 1))
    if (( i % 500 == 0 )); then log "  resolved $i / ${#uniq_accs[@]}"; fi
done

log "Backfilling taxids into merged fasta headers"
FINAL_FAA="$WORK_DIR/kaiju_input.final.faa"
awk -v map="$ACC2TAXID" '
    BEGIN {
        while ((getline line < map) > 0) {
            split(line, a, "\t"); t[a[1]] = a[2]
        }
    }
    /^>/ {
        split(substr($0, 2), parts, "__ACC__")
        n = parts[1]; acc = parts[2]
        tid = (acc in t) ? t[acc] : 1
        printf(">%s_%s\n", n, tid); next
    }
    { gsub(/[^ACDEFGHIKLMNPQRSTVWY]/, "X"); print }
' "$MERGED_FAA" > "$FINAL_FAA"

nseqs=$(grep -c '^>' "$FINAL_FAA" || true)
log "Final merged protein fasta: $nseqs sequences, $(du -h "$FINAL_FAA" | awk '{print $1}')"

# ---- Kaiju BWT + FMI ----------------------------------------------------
BWT_PREFIX="$WORK_DIR/kaiju_dev_db"

log "Building Kaiju BWT"
kaiju_mkbwt -n "$(nproc)" -a ACDEFGHIKLMNPQRSTVWY -o "$BWT_PREFIX" "$FINAL_FAA"

log "Building Kaiju FMI"
kaiju_mkfmi "$BWT_PREFIX"

mv "${BWT_PREFIX}.fmi" "$FINAL_FMI"

log "Done."
du -sh "$DB_DIR" >&2
