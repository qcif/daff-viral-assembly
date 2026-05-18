# Reference Data

## Overview

Reference data is stored in premium blob storage (`daffpremium`) and downloaded to each node's NVMe storage during the pool start task. Two containers are used:

| Container | Who uses it | Contents |
|---|---|---|
| `refdata` | Taxodactyl + wf4 (shared) | `core_nt/` (BLAST), `taxdump/` (taxonkit) |
| `refdata-wf4` | wf4 only | All other wf4 databases |

Both containers are downloaded into `/mnt/nvme/refdata/` on each node, so there is a single mount point. The start task downloads `refdata` first (giving `core_nt/` and `taxdump/`), then copies the *contents* of `refdata-wf4` directly into the same directory, so all databases sit flat alongside each other.

## Reference Databases

The `refdata-wf4` container must be populated before the first run:

| Database | Size | Blob path in `refdata-wf4` | NVMe path |
|---|---|---|---|
| Kraken2 DB | 315GB | `kraken2_db/` | `/mnt/nvme/refdata/kraken2_db/` |
| Kaiju DB | 175GB | `kaiju_db/` | `/mnt/nvme/refdata/kaiju_db/` |
| DIAMOND protein DB | 45GB | `diamond/` | `/mnt/nvme/refdata/diamond/` |
| RVDB taxonomy | 4GB | `rvdb_taxonomy/` | `/mnt/nvme/refdata/rvdb_taxonomy/` |
| GeNomad DB | 1.4GB | `genomad_db/` | `/mnt/nvme/refdata/genomad_db/` |
| Pfam-A.hmm | 1.8GB | `pfam/` | `/mnt/nvme/refdata/pfam/` |
| Contaminant/rRNA refs | 399MB | `rrna_ref/` | `/mnt/nvme/refdata/rrna_ref/` |

Total `refdata-wf4`: ~542GB
Total `refdata` (wf4-relevant): ~270GB
**Grand total on NVMe**: ~812GB (11TB NVMe on L48as_v3 has ample space)

## Uploading Reference Data

### Grant yourself upload permissions

```sh
az login
USER_OBJECT_ID=$(az ad signed-in-user show --query id -o tsv)
SUBSCRIPTION_ID=$(az account show --query id -o tsv)

az role assignment create \
  --role "Storage Blob Data Contributor" \
  --assignee "$USER_OBJECT_ID" \
  --scope "/subscriptions/$SUBSCRIPTION_ID/resourceGroups/daff-biosecurity/providers/Microsoft.Storage/storageAccounts/daffpremium"
```

### Upload with azcopy

Use azcopy for large transfers (much faster than `az storage blob upload-batch`):

```sh
# Upload Kraken2 DB
azcopy copy '/path/to/kraken2_db/*' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/kraken2_db/" \
  --recursive --overwrite=ifSourceNewer

# Upload Kaiju DB
azcopy copy '/path/to/kaiju_db/*' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/kaiju_db/" \
  --recursive --overwrite=ifSourceNewer

# Upload Pfam-A.hmm
azcopy copy '/path/to/Pfam-A.hmm' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/pfam/"

# Upload DIAMOND DB
azcopy copy '/path/to/viral.dmnd' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/diamond/"

# Upload GeNomad DB
azcopy copy '/path/to/genomad_db/*' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/genomad_db/" \
  --recursive --overwrite=ifSourceNewer

# Upload RVDB taxonomy
azcopy copy '/path/to/rvdb_taxonomy/*' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/rvdb_taxonomy/" \
  --recursive --overwrite=ifSourceNewer

# Upload rRNA/contaminant references
azcopy copy '/path/to/rrna_ref/*' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/rrna_ref/" \
  --recursive --overwrite=ifSourceNewer
```

The `core_nt` and `taxdump` data already exist in the `refdata` container from Taxodactyl setup. No re-upload needed.

### Upload test databases for the test pool

For the `view_test` pool, upload small test databases to a `test/` prefix:

```sh
azcopy copy './tests/kraken_test_db/*' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/test/kraken_test_db/" \
  --recursive

azcopy copy './tests/kaiju_test_db/*' \
  "https://daffpremium.blob.core.windows.net/refdata-wf4/test/kaiju_test_db/" \
  --recursive

# These land on the test node at:
#   /mnt/nvme/refdata/test/kraken_test_db/
#   /mnt/nvme/refdata/test/kaiju_test_db/
```

## Generating SAS Tokens

The start task script needs read-access SAS tokens for both containers.

### refdata-wf4 container SAS

With `batch-helpers.sh` (omit blob name for container-level SAS):

```sh
source deploy/azure/batch-helpers.sh
az_sas_generate "" daffpremium refdata-wf4
```

Or directly with Azure CLI:

```sh
set -a && source .env.azure && set +a

EXPIRY=$(date -u -d "+1 year" '+%Y-%m-%dT%H:%MZ')
ACCOUNT_KEY=$(az storage account keys list \
  -g daff-biosecurity -n daffpremium --query '[0].value' -o tsv)

az storage container generate-sas \
  --account-name daffpremium \
  --account-key "$ACCOUNT_KEY" \
  --name refdata-wf4 \
  --permissions rl \
  --expiry "$EXPIRY" \
  --https-only \
  -o tsv
```

### refdata container SAS (shared with taxodactyl)

With `batch-helpers.sh`:

```sh
source deploy/azure/batch-helpers.sh
az_sas_generate "" daffpremium refdata
```

Or directly with Azure CLI:

```sh
az storage container generate-sas \
  --account-name daffpremium \
  --account-key "$ACCOUNT_KEY" \
  --name refdata \
  --permissions rl \
  --expiry "$EXPIRY" \
  --https-only \
  -o tsv
```

Paste both SAS tokens into `deploy/azure/setup.sh.ignore` (gitignored).

## VM Storage Layout (Standard_L48as_v3)

- OS disk: ~30GB at `/`
- Temp disk: ~80GB at `/mnt`
- **NVMe: ~11TB at `/mnt/nvme`** (formatted ext4 by start task)

Available after reference data download: ~10.2TB
