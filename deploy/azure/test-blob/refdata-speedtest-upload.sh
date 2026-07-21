#!/usr/bin/env bash
#
# Upload a test file to both standard and premium blob storage for
# download-speed benchmarking. Uses a throwaway `test-refdata` container
# in each of $STORAGE_ACCOUNT_STD and $STORAGE_ACCOUNT_PREM.
#
# Usage:
#   source .env.azure
#   deploy/azure/refdata-speedtest-upload.sh <path/to/testfile>
#
# Requires: az CLI, azcopy, and .env.azure loaded (STORAGE_ACCOUNT_STD,
# STORAGE_ACCOUNT_PREM, AZURE_STORAGE_ACCOUNT_KEY, NXF_AZURE_REFERENCE_KEY).

set -euo pipefail

source .env.azure

CONTAINER="test-refdata"
SAS_EXPIRY_DAYS=1

src_file="${1:-}"

if [[ -z "$src_file" ]]; then
    echo "Usage: $0 <path/to/testfile>" >&2
    exit 1
fi

if [[ ! -f "$src_file" ]]; then
    echo "ERROR: Source file not found: $src_file" >&2
    exit 1
fi

for var in STORAGE_ACCOUNT_STD STORAGE_ACCOUNT_PREM \
           AZURE_STORAGE_ACCOUNT_KEY NXF_AZURE_REFERENCE_KEY; do
    if [[ -z "${!var:-}" ]]; then
        echo "ERROR: $var is not set. Source .env.azure first." >&2
        exit 1
    fi
done

if ! command -v azcopy >/dev/null 2>&1; then
    echo "ERROR: azcopy is not installed." >&2
    exit 1
fi

blob_name="$(basename "$src_file")"
file_size=$(du -h "$src_file" | cut -f1)

expiry=$(date -u -d "+${SAS_EXPIRY_DAYS} days" '+%Y-%m-%dT%H:%MZ')

echo "Source file:    $src_file ($file_size)"
echo "Blob name:      $blob_name"
echo "Container:      $CONTAINER"
echo "SAS expiry:     $expiry"
echo ""

upload_to_account() {
    local account="$1"
    local account_key="$2"
    local label="$3"

    echo "===== [$label] $account ====="

    echo "Ensuring container '$CONTAINER' exists..."
    az storage container create \
        --account-name "$account" \
        --account-key "$account_key" \
        --name "$CONTAINER" \
        --only-show-errors >/dev/null

    echo "Generating write SAS token..."
    local sas
    sas=$(az storage container generate-sas \
        --account-name "$account" \
        --account-key "$account_key" \
        --name "$CONTAINER" \
        --permissions cw \
        --expiry "$expiry" \
        --https-only \
        --output tsv)

    local url="https://${account}.blob.core.windows.net/${CONTAINER}/${blob_name}?${sas}"

    echo "Uploading with azcopy..."
    azcopy copy "$src_file" "$url"

    echo "[$label] Upload complete: https://${account}.blob.core.windows.net/${CONTAINER}/${blob_name}"
    echo ""
}

upload_to_account "$STORAGE_ACCOUNT_STD"  "$AZURE_STORAGE_ACCOUNT_KEY"  "STANDARD"
upload_to_account "$STORAGE_ACCOUNT_PREM" "$NXF_AZURE_REFERENCE_KEY"    "PREMIUM"

echo "Done. Test file staged in both storage tiers under container '$CONTAINER'."
