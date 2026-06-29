#!/usr/bin/env bash
# =============================================================================
# SLURM job: predict-all inference pass (Phase 3), all 12 combos sequentially.
#
# Submit standalone: sbatch scripts/cluster/run_inference.sh
# Submit after LOO:  sbatch --dependency=afterok:<LOO_JOB_ID> scripts/cluster/run_inference.sh
# (submit_all.sh handles the dependency automatically)
# =============================================================================

#SBATCH --job-name=kinase_inference
#SBATCH --partition=kimq
#SBATCH --cpus-per-task=48
#SBATCH --mem=200G                    # inference holds all 3,228 sites x 420 kinases in RAM
#SBATCH --time=12:00:00
#SBATCH --output=cluster_logs/inference.out
#SBATCH --error=cluster_logs/inference.err

set -euo pipefail

N_JOBS=36
WORKERS=12
DIMENSIONS=64
OUTPUT_DIR=cluster_results

echo "=== Predict-all inference pass ==="
echo "n_jobs=$N_JOBS  workers=$WORKERS  dimensions=$DIMENSIONS"
date

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VENV_DIR="$PROJECT_DIR/cluster_venv"

source "$VENV_DIR/bin/activate"
cd "$PROJECT_DIR"
export PYTHONPATH="$PROJECT_DIR:${PYTHONPATH:-}"

mkdir -p cluster_logs "$OUTPUT_DIR"

python scripts/cluster/cluster_inference_runner.py \
    --dimensions "$DIMENSIONS" \
    --n-jobs     "$N_JOBS"     \
    --workers    "$WORKERS"    \
    --output-dir "$OUTPUT_DIR"

echo "=== Inference done ==="
date
