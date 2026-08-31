#!/usr/bin/env bash
#
# Start task for the `view_speedtest` pool: benchmark azcopy download
# from either standard or premium blob storage.
#
# Before uploading this script to the scripts container, replace the SAS
# token placeholders and un-comment ONE of the two SOURCE_URL lines below
# to select which storage tier to benchmark.

set -x
set -euo pipefail

NVME_MNT="/mnt/nvme"
DEST_DIR="$NVME_MNT/test-refdata"

# Test file name in the `test-refdata` container (matches what was uploaded
# by deploy/azure/refdata-speedtest-upload.sh — set to that file's basename).
TEST_BLOB="k2_core_nt_20251015.tar.gz"

# --- Choose ONE source tier by un-commenting the matching line ---
# Standard blob (daffstandard):
SOURCE_URL="https://daffstandard.blob.core.windows.net/test-refdata/${TEST_BLOB}?se=2026-08-20T00%3A55Z&sp=rl&spr=https&sv=2026-02-06&sr=c&sig=0JIGnkK1x0EMY%2BB%2BR5GmQxsCXOLOcdE0sMeHADFw1HQ%3D"
# Premium blob (daffpremium):
# SOURCE_URL="https://daffpremium.blob.core.windows.net/test-refdata/${TEST_BLOB}?se=2026-08-20T00%3A55Z&sp=rl&spr=https&sv=2026-02-06&sr=c&sig=dfIvQe2d1%2BOVuDyfgG2CIVSvONFwlJdo6C7c81dNwHc%3D"

# -----------------------------------------------------------------------------

detect_nvme_device() {
    local dev
    dev=$(lsblk -ndo NAME,TYPE | awk '$2=="disk" {print $1}' | grep nvme | head -n 1)
    if [[ -z "$dev" ]] && [[ -b /dev/nvme0n1 ]]; then
        dev="nvme0n1"
    fi
    if [[ -z "$dev" ]]; then
        echo "ERROR: No NVMe disk detected." >&2
        exit 1
    fi
    echo "/dev/$dev"
}

mount_nvme() {
    local dev="$1"
    mkdir -p "$NVME_MNT"
    if ! blkid "$dev" >/dev/null 2>&1; then
        mkfs.ext4 -F "$dev"
    fi
    if ! mount | grep -q "$NVME_MNT"; then
        mount "$dev" "$NVME_MNT"
    fi
    df -h "$NVME_MNT" >&2 || true
}

install_azcopy() {
    if command -v azcopy >/dev/null 2>&1; then
        azcopy --version >&2 || true
        return
    fi
    wget -q -O azcopy.tar.gz https://aka.ms/downloadazcopy-v10-linux
    tar -xf azcopy.tar.gz --strip-components=1
    mv azcopy /usr/local/bin/
    chmod +x /usr/local/bin/azcopy
    rm -f azcopy.tar.gz
    azcopy --version >&2 || true
}

format_hms() {
    local s="$1"
    printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}

DOWNLOAD_ELAPSED=0
DOWNLOAD_BYTES=0
DOWNLOAD_MBPS="n/a"
DECOMPRESS_ELAPSED=0
DECOMPRESS_BYTES=0
DECOMPRESS_MBPS="n/a"

benchmark_download() {
    echo "=== Benchmarking azcopy download ===" >&2
    echo "Source: ${SOURCE_URL%%\?*}" >&2
    echo "Dest:   $DEST_DIR" >&2

    mkdir -p "$DEST_DIR"

    local start_ts end_ts
    start_ts=$(date +%s)
    echo "Download start: $(date -u -Iseconds)" >&2

    azcopy copy \
        "${SOURCE_URL}" \
        "${DEST_DIR}/" \
        --check-md5 NoCheck \
        --cap-mbps 0

    end_ts=$(date +%s)
    echo "Download end:   $(date -u -Iseconds)" >&2

    DOWNLOAD_ELAPSED=$((end_ts - start_ts))

    local downloaded_file="$DEST_DIR/$TEST_BLOB"
    if [[ ! -f "$downloaded_file" ]]; then
        echo "ERROR: expected file not found at $downloaded_file" >&2
        ls -lh "$DEST_DIR" >&2 || true
        exit 1
    fi

    DOWNLOAD_BYTES=$(stat -c%s "$downloaded_file")
    if [[ "$DOWNLOAD_ELAPSED" -gt 0 ]]; then
        DOWNLOAD_MBPS=$(awk "BEGIN {printf \"%.2f\", ($DOWNLOAD_BYTES * 8) / ($DOWNLOAD_ELAPSED * 1000000)}")
    fi
}

benchmark_decompress() {
    local downloaded_file="$DEST_DIR/$TEST_BLOB"
    local extract_dir="$DEST_DIR/extracted"

    echo "=== Benchmarking decompression ===" >&2
    echo "Archive: $downloaded_file" >&2
    echo "Extract dir: $extract_dir" >&2

    mkdir -p "$extract_dir"

    local start_ts end_ts
    start_ts=$(date +%s)
    echo "Decompress start: $(date -u -Iseconds)" >&2

    # Handle .tar.gz / .tgz / .tar / .gz — extend as needed.
    case "$TEST_BLOB" in
        *.tar.gz|*.tgz)
            tar -xzf "$downloaded_file" -C "$extract_dir"
            ;;
        *.tar)
            tar -xf "$downloaded_file" -C "$extract_dir"
            ;;
        *.gz)
            gunzip -k -c "$downloaded_file" > "$extract_dir/$(basename "${TEST_BLOB%.gz}")"
            ;;
        *)
            echo "WARN: unknown archive type for $TEST_BLOB — skipping decompression" >&2
            return 0
            ;;
    esac

    end_ts=$(date +%s)
    echo "Decompress end:   $(date -u -Iseconds)" >&2

    DECOMPRESS_ELAPSED=$((end_ts - start_ts))
    DECOMPRESS_BYTES=$(du -sb "$extract_dir" | cut -f1)
    if [[ "$DECOMPRESS_ELAPSED" -gt 0 ]]; then
        DECOMPRESS_MBPS=$(awk "BEGIN {printf \"%.2f\", ($DECOMPRESS_BYTES * 8) / ($DECOMPRESS_ELAPSED * 1000000)}")
    fi

    echo "Extracted contents:" >&2
    du -sh "$extract_dir" >&2 || true
}

print_results() {
    local total=$((DOWNLOAD_ELAPSED + DECOMPRESS_ELAPSED))
    local dl_bytes_h dc_bytes_h
    dl_bytes_h=$(numfmt --to=iec "$DOWNLOAD_BYTES" 2>/dev/null || echo "$DOWNLOAD_BYTES")
    dc_bytes_h=$(numfmt --to=iec "$DECOMPRESS_BYTES" 2>/dev/null || echo "$DECOMPRESS_BYTES")

    echo "" >&2
    echo "===== SPEEDTEST RESULTS =====" >&2
    echo "Tier source:    ${SOURCE_URL%%\?*}" >&2
    echo "Blob:           $TEST_BLOB" >&2
    echo "" >&2
    echo "Download:" >&2
    echo "  Size:         $DOWNLOAD_BYTES bytes ($dl_bytes_h)" >&2
    echo "  Elapsed:      ${DOWNLOAD_ELAPSED}s ($(format_hms $DOWNLOAD_ELAPSED))" >&2
    echo "  Throughput:   ${DOWNLOAD_MBPS} Mbps" >&2
    echo "" >&2
    echo "Decompress:" >&2
    echo "  Extracted:    $DECOMPRESS_BYTES bytes ($dc_bytes_h)" >&2
    echo "  Elapsed:      ${DECOMPRESS_ELAPSED}s ($(format_hms $DECOMPRESS_ELAPSED))" >&2
    echo "  Throughput:   ${DECOMPRESS_MBPS} Mbps (uncompressed out / wall time)" >&2
    echo "" >&2
    echo "Total (download + decompress): ${total}s ($(format_hms $total))" >&2
    echo "=============================" >&2
}

echo "===== START TASK: speedtest =====" >&2
date >&2

NVME_DEV=$(detect_nvme_device)
mount_nvme "$NVME_DEV"
install_azcopy
benchmark_download
benchmark_decompress
print_results

echo "===== START TASK: Completed =====" >&2
exit 0
