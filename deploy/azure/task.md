# Deploying to azure

We have already deployed another workflow, "Taxodactyl" to Azure. Please use its code as a guide for setting up this workflow for Azure:
./taxodactyl

Your task is to create documentation, batch-helpers.sh script and required deployment assets analagous to those in the taxodactyl repo.

In particular, you should create analogs for:

./taxodactyl/.env.azure  # Note some env vars may be different or redundant here
./taxodactyl/.env.azure.sample
./taxodactyl/deployment/azure/.gitignore
./taxodactyl/deployment/azure/pool-setup.json.ignore
./taxodactyl/deployment/azure/run-taxodactyl.sh
./taxodactyl/deployment/azure/pool-setup.json.ignore
./taxodactyl/deployment/azure/pool-setup.json.template
./taxodactyl/deployment/azure/setup.sh.template (the start task)
./taxodactyl/deployment/azure/setup.sh.ignore
./taxodactyl/deployment/azure/storage-policy.json
./taxodactyl/docs/azure/*.md

## Nextflow config

You should also use Taxodactyl's Azure config as inspiration for a Nextflow config for this workflow, though much of it will not apply to this workflow as it will have different env vars and processes:

./taxodactyl/conf/azure.config

It should be written to ./config/azure.config. You probably also want to define ./config/azure.test.config for this workflow, as the full workflow is very cpu/memory demanding.

## Reference data

The following reference data will need to be uploaded to Azure Blob:

contaminant/rRNA refs 399M
Kraken2 DB 315Gb
Kaiju DB 175G
Pfam-A.hmm 1.8G
genomad_db 1.4G
.dmnd DB 45G
rvdb_taxonomy 4G
core_nt 269G
taxonkit 600 MB

Note that core_nt & taxonkit already exist in blob storage for use by Taxodactyl, in a `refdata` bucket. We should upload the additional refdata required here to a new bucket, `refdata-wf4`. The start task for this workflow will need to download both `refdata` and `refdata-wf4` buckets to stage the required data. Use the start task from Taxodactyl as a reference.

## VM type

Initially we will trial production runs on Azure's L48as VM type. That gives us:

- 48 vCPU
- 384GB RAM
- 11TB NVME scratch (for ref data)
- Costs $7.50/hr

However, for development we will use mock databases for Kraken and Kaiju, and try to get away with testing on L4as machines that are much cheaper:

- 4 vCPU
- 32GB RAM
- 960GB NVME scratch
- Costs 65c/hr

## Pool configuration and autoscale

Pool IDs: view and view_test
Autoscale: same as taxodactyl, but scale to 1 node only
Managed identities: I don't think they are required for this workflow.
taskSlotsPerNode: set to cpu count (let NF decide how many tasks to dispatch)
