# Azure Batch Setup — wf4 VIEW Workflow

This directory contains documentation for running the VIEW (Viral Identification and Exploration Workflow) on Azure Batch.

> [!NOTE]
> This workflow shares Azure infrastructure (Batch account, storage accounts) with Taxodactyl. Most of the initial setup only needs to be done once across both workflows. See [01-initial-setup.md](01-initial-setup.md) for what is shared vs workflow-specific.

## Quick Start — running the workflow on Azure

This assumes an Azure Batch pool is already set up (see docs below).

1. Create a `.env.azure` file in the project root (see `.env.azure.sample`)
2. Run the workflow:

```bash
deploy/azure/run-wf4.sh --input index.csv --outdir output/run_001
```

This runs Nextflow locally and submits tasks to an Azure Batch node. With autoscaling enabled, the pool defaults to 0 nodes and provisions an ephemeral node when tasks are submitted. The node is kept warm for 15 minutes after the last task completes so that sequential runs can benefit from re-using the ref data.

**First-run time**: Allow 15-20 minutes for the node to provision and stage ~830GB of reference data from premium blob storage (after the first run, warm nodes skip the download).

---

## Infrastructure Overview

| Resource | Name | Notes |
|---|---|---|
| Batch account | `daffbatch` | Shared with taxodactyl |
| Standard storage | `daffstandard` | Work directory, scripts |
| Premium storage | `daffpremium` | Reference databases |
| Resource group | `daff-biosecurity` | All resources |
| Production pool | `view` | L48as_v3, 48 vCPU, 384GB RAM, 11TB NVMe |
| Test pool | `view_test` | L4as_v3, 4 vCPU, 32GB RAM, 960GB NVMe |

## Reference Data

| Database | Size | Container | Path on NVMe |
|---|---|---|---|
Both blob containers are downloaded into `/mnt/nvme/refdata/` on each node (single mount point).

| Database | Size | Blob container | NVMe path |
|---|---|---|---|
| core_nt (BLAST) | 269GB | `refdata` | `/mnt/nvme/refdata/core_nt/` |
| taxdump (taxonkit) | 600MB | `refdata` | `/mnt/nvme/refdata/taxdump/` |
| Kraken2 DB | 315GB | `refdata-wf4` | `/mnt/nvme/refdata/kraken2_db/` |
| Kaiju DB | 175GB | `refdata-wf4` | `/mnt/nvme/refdata/kaiju_db/` |
| DIAMOND protein DB | 45GB | `refdata-wf4` | `/mnt/nvme/refdata/diamond/` |
| RVDB taxonomy | 4GB | `refdata-wf4` | `/mnt/nvme/refdata/rvdb_taxonomy/` |
| GeNomad DB | 1.4GB | `refdata-wf4` | `/mnt/nvme/refdata/genomad_db/` |
| Pfam-A.hmm | 1.8GB | `refdata-wf4` | `/mnt/nvme/refdata/pfam/` |
| Contaminant/rRNA refs | 399MB | `refdata-wf4` | `/mnt/nvme/refdata/rrna_ref/` |

The `refdata` container is shared with taxodactyl. Only `refdata-wf4` is new.

---

## Documentation Index

1. **[Initial Setup](01-initial-setup.md)** — Azure resources, credentials, quotas
2. **[Reference Data](02-reference-data.md)** — Uploading and staging reference databases
3. **[Start Tasks](03-start-tasks.md)** — Node initialisation and data staging
4. **[Pool Management](04-pool-management.md)** — Creating pools, autoscaling

## Helper Functions

Load the batch helper script for convenient CLI management:

```bash
# Load helpers (automatically loads .env.azure if present)
source deploy/azure/batch-helpers.sh

# View all available commands
az_help

# Common operations
az_pool_list                    # List all pools
az_pool_show                    # Show production pool details
az_pool_show view_test          # Show test pool details
az_pool_resize 1                # Manually scale production pool up to 1 node
az_node_list                    # List nodes in production pool
az_node_logs                    # Download start task logs
az_jobs_list                    # List recent Nextflow jobs
```

## Running the workflow

See [run-wf4.sh](../../deploy/azure/run-wf4.sh)
