#!/bin/bash
#
# Run the wf4 VIEW test case on Azure Batch
#
# Runs the small test case in tests/ (~8MB of mock databases) against the
# view_test pool. Unlike run-wf4.sh this needs NO staged reference data:
# every database is a `path` input, so Nextflow uploads tests/ to the
# az:// workDir and stages it into each task container automatically.
#
# The one exception is the TaxonKit taxdump (~2.6GB), which is staged once per
# node by the pool start task (deploy/azure/setup-wf4-test.sh.template) and
# consumed as a `val` input, so it is never uploaded from this machine.
#
# Override TAXDUMP only if you have staged it somewhere else on the node.
#
# Usage:
#   ./deploy/azure/run-wf4-test.sh [--outdir <dir>] [-resume]
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
OUTDIR="output/test_$RUN_ID"
RESUME=""

# TaxonKit taxdump — a node-local path staged by the pool start task, NOT a
# path on this machine. It is a `val` input, so Nextflow passes it through
# unchanged and the container reads it via the /mnt/nvme/refdata bind mount.
TAXDUMP="${TAXDUMP:-/mnt/nvme/refdata/taxdump}"

# Parse command line arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
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
            echo "Usage: $0 [--outdir <dir>] [-resume]"
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
echo -e "${YELLOW}=== VIEW Test Case Azure Batch Configuration ===${NC}"
echo "Input:           tests/index_test.csv (from test profile)"
echo "Output dir:      $OUTDIR"
echo "Test databases:  tests/ (staged to the node by Nextflow)"
echo "Taxdump:         $TAXDUMP (on Azure Batch nodes)"
echo "Pool:            view_test"
echo "Profile:         azure_test,test"
echo "Resume:          ${RESUME:-false}"
echo ""

# Confirm execution
read -p "Continue with test execution? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    echo "Execution cancelled"
    exit 0
fi

echo ""
echo -e "${GREEN}=== Starting VIEW Test Case ===${NC}"
echo ""

mkdir -p "$OUTDIR"

# The trailing `test` profile supplies the tests/ database paths AND is what
# triggers GENOMAD_DOWNLOAD_DB — see the note in conf/azure.test.config.
nextflow run main.nf \
    -profile azure_test,test \
    --outdir "$OUTDIR" \
    --taxdump "$TAXDUMP" \
    --analyst_name "${ANALYST_NAME:-Tester}" \
    --facility "${FACILITY_NAME:-Unknown}" \
    $RESUME

exit_code=$?

echo ""
if [[ $exit_code -eq 0 ]]; then
    echo -e "${GREEN}=== Test Case Completed Successfully ===${NC}"
    echo ""
    echo "Output directory: $OUTDIR"
else
    echo -e "${RED}=== Test Case Failed ===${NC}"
    echo ""
    echo "Exit code: $exit_code"
    echo "Check .nextflow.log for details"
fi

exit $exit_code
