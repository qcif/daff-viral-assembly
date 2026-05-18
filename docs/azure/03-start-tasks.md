# Start Tasks

## Overview

The start task runs automatically when a pool node is provisioned. For wf4 it:

1. Detects and formats the NVMe device (`/dev/nvme0n1`)
2. Mounts NVMe to `/mnt/nvme`
3. Downloads and installs azcopy
4. Downloads the `refdata` container → `/mnt/nvme/refdata` (shared with taxodactyl: core_nt, taxdump)
5. Downloads the `refdata-wf4` container → `/mnt/nvme/refdata-wf4` (wf4-specific databases)

Warm-node detection: if the target directories already exist, downloads are skipped. This is important for the autoscale "keep warm" window — nodes that don't deallocate between jobs skip the ~25 min download.

## Script Files

| File | Purpose |
|---|---|
| `deploy/azure/setup.sh.template` | Version-controlled template with `<PLACEHOLDER>` tokens |
| `deploy/azure/setup.sh.ignore` | Gitignored copy with real SAS tokens (used for deployment) |

## Populating SAS Tokens

1. Generate SAS tokens (see [02-reference-data.md](02-reference-data.md))
2. Copy the template: `cp deploy/azure/setup.sh.template deploy/azure/setup.sh.ignore`
3. Replace `<REFDATA_SAS_TOKEN>` and `<REFDATA_WF4_SAS_TOKEN>` in `setup.sh.ignore`

The BLOB_URLs in the script should look like:

```bash
REFDATA_BLOB_URL="https://daffpremium.blob.core.windows.net/refdata?se=2027-...&sp=rl&..."
REFDATA_WF4_BLOB_URL="https://daffpremium.blob.core.windows.net/refdata-wf4?se=2027-...&sp=rl&..."
```

## Deploying the Start Task Script

Upload `setup.sh.ignore` to blob storage so Batch nodes can download it at startup:

```sh
source deploy/azure/batch-helpers.sh
az_storage_upload deploy/azure/setup.sh.ignore setup-wf4.sh
```

Generate a SAS token for the uploaded script:

```sh
az_sas_generate setup-wf4.sh
```

Copy the full URL (with SAS token) into `deploy/azure/pool-setup.json.ignore`:

```json
"startTask": {
  "resourceFiles": [{
    "httpUrl": "https://daffstandard.blob.core.windows.net/scripts/setup-wf4.sh?<SAS>",
    "filePath": "setup.sh"
  }]
}
```

## Applying the Configuration

```sh
# Update pool start task
az_pool_update --json deploy/azure/pool-setup.json.ignore view

# Force existing nodes to re-run the start task (recreate)
az_pool_resize 0 view --yes && az_pool_resize 1 view
```

## Debugging Start Task Failures

```sh
source deploy/azure/batch-helpers.sh
az_load_env

# Check node state and start task result
az_node_list

# Download and view stderr (verbose due to set -x)
az_node_logs               # production pool stderr
az_node_logs view_test     # test pool stderr
az_node_logs view stdout   # stdout (mkfs output)
```

Common failure causes:

| Symptom | Likely cause |
|---|---|
| `azcopy: command not found` | wget failed — check internet access / proxy |
| `403` / `AuthenticationFailed` | SAS token expired or incorrect permissions |
| `No such file or directory` on refdata path | Download completed to wrong directory; check azcopy output |
| Node stuck in `starting` | Start task still running (first-run download takes 20–35 min) |

## Performance

| Step | Duration |
|---|---|
| NVMe format + mount | ~10 seconds |
| azcopy install | ~15 seconds |
| Download `refdata` (~270GB) | ~3 minutes |
| Download `refdata-wf4` (~542GB) | ~6 minutes |
| **Total (cold node)** | **~10 minutes** |

Throughput: ~1.5 GB/s (premium blob → NVMe via azcopy parallel transfer)
