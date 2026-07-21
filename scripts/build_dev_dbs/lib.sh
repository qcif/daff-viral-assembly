#!/usr/bin/env bash
# Shared helpers for the build_dev_dbs scripts.

set -euo pipefail

# Docker images for kraken2 / kaiju. Override via env if needed.
KRAKEN2_IMAGE="${KRAKEN2_IMAGE:-staphb/kraken2:2.1.3}"
KAIJU_IMAGE="${KAIJU_IMAGE:-nanozoo/kaiju:1.10.1--90efeeb}"

log() { printf '[%s] %s\n' "$(date +%H:%M:%S)" "$*" >&2; }

require_cmd() {
    for cmd in "$@"; do
        command -v "$cmd" >/dev/null 2>&1 || {
            echo "ERROR: required command not found: $cmd" >&2
            exit 1
        }
    done
}

# Ensure a docker image is present locally (pull once).
ensure_image() {
    local img="$1"
    if ! docker image inspect "$img" >/dev/null 2>&1; then
        log "Pulling docker image: $img"
        docker pull "$img" >&2
    fi
}

# Run a command inside a docker image, mounting $OUT_ROOT so all paths under
# it resolve identically inside and outside the container. Runs as the current
# uid/gid so outputs are owned by the invoking user.
# Usage: docker_run <image> <cmd> [args...]
docker_run() {
    local img="$1"; shift
    docker run --rm \
        -u "$(id -u):$(id -g)" \
        -v "$OUT_ROOT:$OUT_ROOT" \
        -w "$OUT_ROOT" \
        "$img" "$@"
}

kraken2_build() { docker_run "$KRAKEN2_IMAGE" kraken2-build "$@"; }
kaiju_mkbwt()   { docker_run "$KAIJU_IMAGE" kaiju-mkbwt "$@"; }
kaiju_mkfmi()   { docker_run "$KAIJU_IMAGE" kaiju-mkfmi "$@"; }

# Parse seed_taxa.txt into per-section accession lists.
# Usage: parse_section <section> <seed_file>
parse_section() {
    local section="$1" file="$2"
    awk -v want="[$section]" '
        /^\[.*\]$/ { in_section = ($0 == want); next }
        in_section && /^[^#[:space:]]/ { print $1 }
    ' "$file"
}

# efetch a single nucleotide accession to a fasta file (idempotent).
efetch_nuccore() {
    local acc="$1" out="$2"
    if [[ -s "$out" ]]; then
        log "  cached: $out"
        return 0
    fi
    log "  fetching nuccore: $acc"
    curl -fsSL \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
        -o "$out.tmp"
    mv "$out.tmp" "$out"
}

# efetch protein sequences linked to a nucleotide accession (idempotent).
efetch_linked_proteins() {
    local acc="$1" out="$2"
    if [[ -s "$out" ]]; then
        log "  cached: $out"
        return 0
    fi
    log "  fetching linked proteins: $acc"
    local ids
    ids=$(curl -fsSL \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/elink.fcgi?dbfrom=nuccore&db=protein&id=${acc}&idtype=acc" \
        | grep -oE '<Id>[0-9]+</Id>' \
        | sed -E 's|</?Id>||g' \
        | paste -sd, -)
    if [[ -z "$ids" ]]; then
        log "  (no linked proteins for $acc)"
        : > "$out"
        return 0
    fi
    curl -fsSL \
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=${ids}&rettype=fasta&retmode=text" \
        -o "$out.tmp"
    mv "$out.tmp" "$out"
}

# Resolve OUT_ROOT. Precedence:
#   1. $OUT_ROOT env var
#   2. $REPO_ROOT/refdata_dev if running inside the wf4 repo checkout
resolve_out_root() {
    if [[ -n "${OUT_ROOT:-}" ]]; then
        mkdir -p "$OUT_ROOT"
        OUT_ROOT="$(cd -- "$OUT_ROOT" && pwd)"
        return 0
    fi
    local here="$1"
    local repo_root
    repo_root="$(cd -- "$here/../.." && pwd)"
    if [[ -f "$repo_root/nextflow.config" ]]; then
        OUT_ROOT="$repo_root/refdata_dev"
        mkdir -p "$OUT_ROOT"
        return 0
    fi
    echo "ERROR: OUT_ROOT not set and not running inside wf4 repo checkout." >&2
    echo "       Set OUT_ROOT=/path/to/writable/dir and rerun." >&2
    exit 1
}
