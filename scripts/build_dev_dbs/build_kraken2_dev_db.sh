#!/usr/bin/env bash
# Build a small "dev" Kraken2 database (~500 MB target) for wf4.
# Idempotent: rerun to resume after failures. Delete $OUT_ROOT/kraken2_db to
# force a full rebuild.
#
# Env overrides:
#   OUT_ROOT         Output root (default: <repo>/refdata_dev). Must be writable.
#   KRAKEN2_IMAGE    Docker image (default: staphb/kraken2:2.1.3)

set -euo pipefail

HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
# shellcheck source=lib.sh
source "$HERE/lib.sh"

resolve_out_root "$HERE"

SEED_FILE="$HERE/seed_taxa.txt"
DB_DIR="$OUT_ROOT/kraken2_db"
WORK_DIR="$OUT_ROOT/.build/kraken2"
FASTA_DIR="$WORK_DIR/fastas"

SIZE_BUDGET_MB=500

require_cmd docker curl awk
ensure_image "$KRAKEN2_IMAGE"

log "OUT_ROOT = $OUT_ROOT"

if [[ -f "$DB_DIR/hash.k2d" && -f "$DB_DIR/opts.k2d" && -f "$DB_DIR/taxo.k2d" ]]; then
    log "Kraken2 dev DB already present at $DB_DIR — nothing to do."
    du -sh "$DB_DIR" >&2
    exit 0
fi

mkdir -p "$DB_DIR" "$FASTA_DIR"

log "Downloading NCBI taxonomy into $DB_DIR/taxonomy (skips if present)"
kraken2_build --download-taxonomy --db "$DB_DIR"

log "Downloading viral RefSeq library"
kraken2_build --download-library viral --db "$DB_DIR"

log "Fetching per-accession seed sequences (fabaceae + context)"
for section in fabaceae context; do
    while IFS= read -r acc; do
        [[ -z "$acc" ]] && continue
        out="$FASTA_DIR/${acc}.fa"
        efetch_nuccore "$acc" "$out"
        log "  adding to library: $acc"
        kraken2_build --add-to-library "$out" --db "$DB_DIR" --no-masking
    done < <(parse_section "$section" "$SEED_FILE")
done

lib_bytes=$(du -sb "$DB_DIR/library" | awk '{print $1}')
lib_mb=$(( lib_bytes / 1024 / 1024 ))
log "Library staged: ${lib_mb} MB (budget ${SIZE_BUDGET_MB} MB for final DB)"
if (( lib_mb > SIZE_BUDGET_MB * 4 )); then
    log "WARN: library is >4x budget — final DB will likely exceed ${SIZE_BUDGET_MB} MB."
    log "      Trim seed_taxa.txt or lower --kmer-len."
fi

log "Building Kraken2 index (this is the slow step)"
kraken2_build --build --db "$DB_DIR" --threads "$(nproc)"

log "Cleaning intermediate library files"
kraken2_build --clean --db "$DB_DIR"

log "Done."
du -sh "$DB_DIR" >&2
