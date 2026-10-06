#!/bin/bash
#
# Run the MEGABLAST batch-vs-per-query performance test on Azure Batch.
#
# Each workflow is its own invocation, so the `view` pool can scale back
# to 0 nodes between them rather than sharing one warm node. Each run
# provisions its own node (the NODE_UP task inside batch.nf/per_query.nf),
# runs its BLAST workload, and publishes under results/<run_id>/<workflow>/.
#
# This uses the production `view` pool. Check manually that it's not busy
# before running this.
#
# Usage:
#   ./tests/blast-perf/run.sh batch                  # starts a new run_id, prints it
#   ./tests/blast-perf/run.sh per_query <run_id>      # reuses that run_id
#   ./tests/blast-perf/run.sh report <run_id>         # generates REPORT.md from both
#
# Run from the repository root. Wait for the pool to scale back to 0
# nodes (about 15 minutes after the queue empties) before starting the
# second workflow, so it doesn't inherit the first run's warm node/cache.

set -euo pipefail

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m'

TEST_DIR="tests/blast-perf"
CONFIG="${TEST_DIR}/blast-perf.config"

usage() {
    echo "Usage: $0 <batch|per_query> [run_id]"
    echo "       $0 report <run_id>"
    exit 1
}

MODE="${1:-}"
[[ -z "$MODE" ]] && usage

if [[ ! -f main.nf ]]; then
    echo -e "${RED}ERROR: run from the repository root (main.nf not found)${NC}"
    exit 1
fi

if [[ "$MODE" == "report" ]]; then
    RUN_ID="${2:-}"
    [[ -z "$RUN_ID" ]] && usage
    RESULTS_DIR="${TEST_DIR}/results/${RUN_ID}"

    (
        cd "$TEST_DIR"
        ../../venv/bin/python3 summarise.py \
            "results/${RUN_ID}/batch" \
            "results/${RUN_ID}/per_query"
        mv REPORT.md "results/${RUN_ID}/REPORT.md"
    )
    echo "Report: ${RESULTS_DIR}/REPORT.md"
    exit 0
fi

if [[ "$MODE" != "batch" && "$MODE" != "per_query" ]]; then
    usage
fi

RUN_ID="${2:-$(date +"%Y%m%d_%H%M%S")}"
RESULTS_DIR="${TEST_DIR}/results/${RUN_ID}/${MODE}"

if [[ ! -f .env.azure ]]; then
    echo -e "${RED}ERROR: .env.azure not found${NC}"
    exit 1
fi

echo -e "${YELLOW}Loading Azure credentials from .env.azure${NC}"
set -a
source .env.azure
set +a
source deploy/azure/batch-helpers.sh

echo ""
echo -e "${YELLOW}=== BLAST batch-vs-per-query performance test: ${MODE} ===${NC}"
echo "Run ID:      $RUN_ID"
echo "Results dir: $RESULTS_DIR"
echo "Pool:        view"
echo ""

read -p "This runs real work on the production 'view' pool and costs money (~\$4-6 for this one workflow). Make sure it's not already busy. Continue? (yes/no): " confirm
if [[ "$confirm" != "yes" ]]; then
    echo "Execution cancelled"
    exit 0
fi

mkdir -p "$RESULTS_DIR"

echo ""
echo -e "${GREEN}=== Running ${MODE}.nf ===${NC}"
nextflow run "${TEST_DIR}/${MODE}.nf" \
    -c "$CONFIG" \
    --outdir "$RESULTS_DIR"

echo ""
echo -e "${GREEN}=== Done: ${MODE} ===${NC}"
echo "Run ID: $RUN_ID"
echo ""
echo "The 'view' pool will scale back to 0 nodes about 15 minutes after" \
    "this job queue empties."
echo ""
if [[ "$MODE" == "batch" ]]; then
    echo "Wait for the pool to reach 0 nodes, then run:"
    echo "  ./tests/blast-perf/run.sh per_query $RUN_ID"
else
    echo "Once both have run, generate the report with:"
    echo "  ./tests/blast-perf/run.sh report $RUN_ID"
fi
