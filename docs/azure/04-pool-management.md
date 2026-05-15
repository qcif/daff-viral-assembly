# Pool Management

## Pool Overview

Two pools are used for this workflow:

| Pool | VM SKU | vCPU | RAM | NVMe | Cost | Use |
|---|---|---|---|---|---|---|
| `view` | Standard_L48as_v3 | 48 | 384GB | 11TB | $7.50/hr | Production |
| `view_test` | Standard_L4as_v3 | 4 | 32GB | 960GB | $0.65/hr | Development / testing |

> [!IMPORTANT]
> Do not use Nextflow's `autoPoolMode`. The pool must be created manually so that the start task can stage reference data. Auto-created pools have no start task and no reference data.

Managed identities are not required for this workflow (no Key Vault).

## Creating the Production Pool

```sh
source deploy/azure/batch-helpers.sh  # loads .env.azure automatically

# Create pool with autoscaling
az_pool_create --json deploy/azure/pool-setup.json.ignore --autoscale
```

Or with the raw Azure CLI:

```sh
set -a && source .env.azure && set +a

az batch pool create \
  --account-name $AZURE_BATCH_ACCOUNT_NAME \
  --account-endpoint $AZURE_BATCH_ENDPOINT \
  --json-file deploy/azure/pool-setup.json.ignore
```

## Creating the Test Pool

```sh
az_pool_create --json deploy/azure/pool-setup-test.json.ignore --autoscale
```

## Autoscale Formula

Both pools use the same autoscale formula (defined in `batch-helpers.sh`):

- Scale to **1 node** if any tasks are pending in the last 15 minutes
- Scale to **0 nodes** when idle for 15 minutes
- Evaluation interval: 5 minutes (minimum allowed by Azure)

```
$keepWarmMinutes = 15;
$concurrentNodes = 1;
interval = TimeInterval_Minute * $keepWarmMinutes;
$tasks = max($PendingTasks.GetSample(1), avg($PendingTasks.GetSample(interval)));
$TargetDedicatedNodes = $tasks > 0 ? $concurrentNodes : 0;
$NodeDeallocationOption = taskcompletion;
```

## Enabling Autoscale Manually

```sh
# Enable autoscale on production pool
az_pool_update --autoscale view

# Enable autoscale on test pool
az_pool_update --autoscale view_test
```

## Pool Management Commands

```sh
# List all pools
az_pool_list

# Show pool details
az_pool_show           # production
az_pool_show view_test # test

# Resize pools
az_pool_resize 1                # scale production to 1 node
az_pool_resize 0                # scale production to 0
az_pool_resize 1 view_test      # scale test to 1 node

# Delete a pool (with confirmation)
az_pool_delete view

# Update start task (recreate node to apply)
az_pool_update --json deploy/azure/pool-setup.json.ignore view
az_pool_resize 0 --yes && az_pool_resize 1
```

## taskSlotsPerNode

Both pools set `taskSlotsPerNode` equal to the node's vCPU count. This allows Nextflow to dispatch as many concurrent tasks as the node has CPUs, letting Nextflow control scheduling via `maxForks` in `conf/azure.config`.

- `view`: `taskSlotsPerNode=48`
- `view_test`: `taskSlotsPerNode=4`

## Node Timing

**Cold start (new node)**: ~30–40 minutes
- ~5 min: VM provisioning
- ~25–35 min: start task (NVMe format + mount + download ~830GB at ~1.5 GB/s)

**Warm node (data already staged)**: ~3–5 minutes (just VM provisioning + start task skips download)

## Next Steps

- [Back to Initial Setup](01-initial-setup.md)
- [Reference Data](02-reference-data.md)
- [Start Tasks](03-start-tasks.md)
