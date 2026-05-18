# Initial Setup

> [!NOTE]
> The Azure infrastructure (Batch account, storage accounts, resource group) is shared with Taxodactyl and only needs to be created once. If it already exists, skip to **Step 3: Get Credentials**.

## Prerequisites

- Azure CLI installed and configured
- Access to the `DAFF Biosecurity` subscription
- `source deploy/azure/batch-helpers.sh` loaded for helper functions

## Step 1: Azure Resources (shared, one-time setup)

These resources already exist if Taxodactyl is deployed. Skip if present.

```sh
SUBSCRIPTION="DAFF Biosecurity"
REGION=australiaeast
RESOURCE_GROUP=daff-biosecurity
STORAGE_ACCOUNT_PREM=daffpremium
STORAGE_ACCOUNT_STD=daffstandard
BATCH_ACCOUNT=daffbatch

az account set --subscription "$SUBSCRIPTION"
az group create -n $RESOURCE_GROUP -l $REGION

# Premium storage for reference data
az storage account create \
  -n $STORAGE_ACCOUNT_PREM -g $RESOURCE_GROUP -l $REGION \
  --sku Premium_LRS --kind BlockBlobStorage

# Standard storage for work directory and scripts
az storage account create \
  -n $STORAGE_ACCOUNT_STD -g $RESOURCE_GROUP -l $REGION \
  --sku Standard_LRS

# Batch account
az batch account create \
  --name $BATCH_ACCOUNT \
  --resource-group $RESOURCE_GROUP \
  --location $REGION \
  --storage-account $STORAGE_ACCOUNT_STD
```

## Step 2: Create wf4-specific Storage Container

The `refdata-wf4` container in premium storage is new (not shared with Taxodactyl):

```sh
az storage container create \
  --name refdata-wf4 \
  --account-name $STORAGE_ACCOUNT_PREM
```

Ensure the work and scripts containers exist (shared with Taxodactyl):

```sh
az storage container create --name workdata --account-name $STORAGE_ACCOUNT_STD
az storage container create --name scripts  --account-name $STORAGE_ACCOUNT_STD
```

Apply the storage lifecycle policy to auto-clean work files after 14 days:

```sh
az storage account management-policy create \
  --account-name $STORAGE_ACCOUNT_STD \
  --resource-group $RESOURCE_GROUP \
  --policy deploy/azure/storage-policy.json
```

## Step 3: Get Credentials

```sh
# Retrieve keys (store in 1Password and add to .env.azure)
az storage account keys list -g daff-biosecurity -n daffstandard
az storage account keys list -g daff-biosecurity -n daffpremium
az batch account keys list -g daff-biosecurity -n daffbatch
```

Add to `.env.azure` (copy from `.env.azure.sample`):

```bash
AZURE_BATCH_ACCESS_KEY="<batch-primary-key>"
AZURE_STORAGE_ACCOUNT_KEY="<standard-storage-key1>"
NXF_AZURE_REFERENCE_KEY="<premium-storage-key1>"
```

## Step 4: Request VM Quota

The production pool uses `Standard_L48as_v3`. Request quota if not already available:

1. Sign into Azure Portal → search "Quotas"
2. Request quota for "Batch" type
3. Select "LS series" and "Pools per Batch account"
4. Request at least 1 pool and 1 `L48as_v3` node

Check available SKUs in your region:

```sh
az batch location list-skus --location australiaeast | grep -i l48
```

## Next Steps

1. [Upload reference data](02-reference-data.md)
2. [Configure start tasks](03-start-tasks.md)
3. [Create batch pools](04-pool-management.md)
