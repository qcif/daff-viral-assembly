#!/bin/bash
#
# Run the wf4 VIEW workflow on Azure Batch
#
# This script runs the VIEW Nextflow workflow using Azure Batch for
# compute-intensive processes (Kraken2, Kaiju, BLAST, SPAdes, etc.).
# Reference data is pre-staged to NVMe storage on Azure Batch nodes.
#
# Usage:
#   ./deploy/azure/run-wf4.sh \
#       --params-file params/azure_params.yml
#       --input /path/to/index.csv \
#       --outdir /path/to/output
#

set -euo pipefail

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PID=$$
RUN_ID="$(date +"%Y%m%d_%H%M%S")_$PID"

# Default values
INPUT="index.csv"
PARAMS_FILE="params/azure_params.yml"
OUTDIR="output/$RUN_ID"
RESUME="false"

# Azure Batch node paths (staged by start task to NVMe under /mnt/nvme/refdata/)
# Shared refdata (also used by taxodactyl):
# BLASTN_DB="/mnt/nvme/refdata/core_nt/core_nt"
# TAXDUMP="/mnt/nvme/refdata/taxdump/taxdump"
# # wf4-specific refdata (merged directly into /mnt/nvme/refdata/):
# KRAKEN2_DB="/mnt/nvme/refdata/kraken2_db"
# # Path to the .fmi itself — the workflow derives the containing directory and
# # globs it for *.fmi, *names.dmp and *nodes.dmp
# KAIJU_DB="/mnt/nvme/refdata/kaiju_db/kaiju_db.fmi"
# HMMER_DB="/mnt/nvme/refdata/pfam/Pfam-A.hmm"
# PROT_DB="/mnt/nvme/refdata/diamond/viral.dmnd"
# GENOMAD_DB="/mnt/nvme/refdata/genomad_db"
# RVDB_TAXONOMY="/mnt/nvme/refdata/rvdb_taxonomy"
# RRNA_REF="/mnt/nvme/refdata/rrna_ref"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --params-file)
            if [[ $# -lt 2 ]]; then
                echo -e "${RED}ERROR: --params-file requires a value${NC}"
                exit 1
            fi

            PARAMS_FILE="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        --input)
            INPUT="$2"
            shift 2
            ;;
        -resume)
            RESUME=true
            shift
            ;;
        -h|--help)
            echo "Usage:"
            echo "  $0 [--params-file <params.yml>] [--input <index.csv>] [--outdir <outdir>] [-resume]"
            echo ""
            echo "Requirements:"
            echo "  --params-file FILE   Nextflow YAML parameters file"
            echo "                       Default: params/azure_params.yml"
            echo "Options:"
            echo "  --input FILE         Input index CSV file"
            echo "                       Default: index.csv"
            echo "  --outdir DIR         Output directory"
            echo "                       Default: output/$RUN_ID"
            echo "  -resume              Resume the previous Nextflow execution"
            echo "                       Default: false"
            echo "  -h, --help           Show this help message"
            exit 0
            ;;

        *)
            echo -e "${RED}ERROR: Unknown argument: $1${NC}"
            echo ""
            echo "Usage: $0 [--params-file <params.yml>] [--input <index.csv>] [--outdir <outdir>] [-resume]"
            exit 1
            ;;
    esac
done


# ---------------------------------------------------------------------------
# Validate required files
# ---------------------------------------------------------------------------

if [[ ! -f "$PARAMS_FILE" ]]; then
    echo -e "${RED}ERROR: Params file not found: $PARAMS_FILE${NC}"
    exit 1
fi


if [[ ! -f "main.nf" ]]; then
    echo -e "${RED}ERROR: main.nf not found${NC}"
    echo "Run this script from the repository root directory."
    exit 1
fi

# Load Azure credentials
if [[ ! -f .env.azure ]]; then
    echo -e "${RED}ERROR: .env.azure not found${NC}"
    echo "Please run this script from the repository root directory"
    exit 1
fi

echo -e "${YELLOW}Loading Azure credentials from .env.azure${NC}"
set -a
source .env.azure
set +a

# Verify credentials are set
if [[ -z "${AZURE_STORAGE_ACCOUNT_KEY:-}" ]]; then
    echo -e "${RED}ERROR: AZURE_STORAGE_ACCOUNT_KEY not set in .env.azure${NC}"
    exit 1
fi

if [[ -z "${AZURE_BATCH_ACCESS_KEY:-}" ]]; then
    echo -e "${RED}ERROR: AZURE_BATCH_ACCESS_KEY not set in .env.azure${NC}"
    exit 1
fi

# Show configuration
echo ""
echo -e "${YELLOW}=== VIEW Workflow Azure Batch Run ===${NC}"
echo "Profile:       azure"
echo "Params file:   $PARAMS_FILE"
echo "Resume:        $RESUME"
echo ""
echo "Workflow parameters will be read from:"
echo "  $PARAMS_FILE"
echo ""

# Confirm execution
read -p "Continue with workflow execution? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    echo "Execution cancelled"
    exit 0
fi

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# Construct Nextflow command
# ---------------------------------------------------------------------------

nextflow_args=(
    run
    main.nf
    -profile
    azure
    -params-file
    "$PARAMS_FILE"
)

if [[ "$RESUME" == true ]]; then
    nextflow_args+=("-resume")
fi
    
# ---------------------------------------------------------------------------
# Run workflow
# ---------------------------------------------------------------------------

echo ""
echo -e "${GREEN}=== Starting VIEW Workflow ===${NC}"
echo ""

set +e
nextflow "${nextflow_args[@]}"
exit_code=$?
set -e

# ---------------------------------------------------------------------------
# Report result
# ---------------------------------------------------------------------------


echo ""
if [[ $exit_code -eq 0 ]]; then
    echo -e "${GREEN}=== Workflow Completed Successfully ===${NC}"
    echo ""
    echo "Output directory: $OUTDIR"
else
    echo -e "${RED}=== Workflow Failed ===${NC}"
    echo ""
    echo "Exit code: $exit_code"
    echo "Check .nextflow.log for details"
fi

exit $exit_code
