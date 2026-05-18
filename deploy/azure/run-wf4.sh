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
OUTDIR="output/$RUN_ID"
RESUME=""

# Azure Batch node paths (staged by start task to NVMe under /mnt/nvme/refdata/)
# Shared refdata (also used by taxodactyl):
BLASTN_DB="/mnt/nvme/refdata/core_nt/core_nt"
TAXDUMP="/mnt/nvme/refdata/taxdump"
# wf4-specific refdata (merged directly into /mnt/nvme/refdata/):
KRAKEN2_DB="/mnt/nvme/refdata/kraken2_db"
KAIJU_DB_PATH="/mnt/nvme/refdata/kaiju_db"
KAIJU_DBNAME="kaiju_db.fmi"
HMMER_DB="/mnt/nvme/refdata/pfam/Pfam-A.hmm"
PROT_DB="/mnt/nvme/refdata/diamond/viral.dmnd"
GENOMAD_DB="/mnt/nvme/refdata/genomad_db"
RVDB_TAXONOMY="/mnt/nvme/refdata/rvdb_taxonomy"
RRNA_REF="/mnt/nvme/refdata/rrna_ref"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --input)
            INPUT="$2"
            shift 2
            ;;
        --outdir)
            OUTDIR="$2"
            shift 2
            ;;
        -resume)
            RESUME="-resume"
            shift
            ;;
        *)
            echo -e "${RED}ERROR: Unknown argument: $1${NC}"
            echo ""
            echo "Usage: $0 --input <index.csv> --outdir <dir> [-resume]"
            exit 1
            ;;
    esac
done

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
echo -e "${YELLOW}=== VIEW Workflow Azure Batch Configuration ===${NC}"
echo "Input:           $INPUT"
echo "Output dir:      $OUTDIR"
echo "BLAST DB:        $BLASTN_DB (on Azure Batch nodes)"
echo "Taxdump:         $TAXDUMP (on Azure Batch nodes)"
echo "Kraken2 DB:      $KRAKEN2_DB (on Azure Batch nodes)"
echo "Kaiju DB:        $KAIJU_DB_PATH/$KAIJU_DBNAME (on Azure Batch nodes)"
echo "Pfam HMM:        $HMMER_DB (on Azure Batch nodes)"
echo "Protein DB:      $PROT_DB (on Azure Batch nodes)"
echo "GeNomad DB:      $GENOMAD_DB (on Azure Batch nodes)"
echo "RVDB taxonomy:   $RVDB_TAXONOMY (on Azure Batch nodes)"
echo "rRNA ref:        $RRNA_REF (on Azure Batch nodes)"
echo "Profile:         azure"
echo "Resume:          ${RESUME:-false}"
echo ""

# Confirm execution
read -p "Continue with workflow execution? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    echo "Execution cancelled"
    exit 0
fi

echo ""
echo -e "${GREEN}=== Starting VIEW Workflow ===${NC}"
echo ""

mkdir -p "$OUTDIR"

# Run Nextflow with Azure Batch profile
nextflow run main.nf \
    -profile azure \
    --input "$INPUT" \
    --outdir "$OUTDIR" \
    --blastn_db "$BLASTN_DB" \
    --taxdump "$TAXDUMP" \
    --kraken2_db "$KRAKEN2_DB" \
    --kaiju_db_path "$KAIJU_DB_PATH" \
    --kaiju_dbname "$KAIJU_DBNAME" \
    --hmmer_db "$HMMER_DB" \
    --prot_db "$PROT_DB" \
    --genomad_db "$GENOMAD_DB" \
    --rvdb_taxonomy "$RVDB_TAXONOMY" \
    --rrna_ref "$RRNA_REF" \
    --analyst_name "${ANALYST_NAME:-}" \
    --facility "${FACILITY_NAME:-}" \
    $RESUME

exit_code=$?

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
