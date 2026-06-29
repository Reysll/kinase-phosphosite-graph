#!/usr/bin/env bash
# =============================================================================
# Master submission script: freeze trials, then submit all SLURM jobs in order.
#
# BEFORE RUNNING THIS:
#   1. Transfer the project to the cluster (rsync or git clone)
#   2. Run setup once: bash scripts/cluster/setup_env.sh
#   3. Edit run_loo.sh, run_inference.sh, run_analysis.sh:
#        --partition=<your partition>
#        --account=<your account>
#   4. Run this script from the project root on a login node:
#        bash scripts/cluster/submit_all.sh
#
# The script submits 3 job stages with SLURM dependencies so they run
# in order automatically.
#
# Outputs land in:
#   cluster_results/{experiment_name}/     <- LOO results (12 combos)
#   cluster_results/inference_all/         <- predict-all inference (12 combos)
#   outputs_prediction/analysis/           <- figures + heatmaps
#   cluster_logs/                          <- SLURM stdout/stderr
# =============================================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
cd "$PROJECT_DIR"

VENV_DIR="$PROJECT_DIR/cluster_venv"
if [ ! -f "$VENV_DIR/bin/activate" ]; then
    echo "ERROR: cluster_venv not found. Run setup_env.sh first."
    exit 1
fi

mkdir -p cluster_logs cluster_results

# --------------------------------------------------------------------------
# Step 0: freeze the full all-sites trial set (runs interactively — fast)
# --------------------------------------------------------------------------
echo "=== Step 0: Freezing full trial set (all sites) ==="
source "$VENV_DIR/bin/activate"
export PYTHONPATH="$PROJECT_DIR:${PYTHONPATH:-}"

if [ -f "outputs_prediction/frozen_trials_all.csv.gz" ]; then
    echo "frozen_trials_all.csv.gz already exists — skipping freeze."
else
    python -m src_prediction.runners.run_freeze_all_trials
fi

deactivate 2>/dev/null || true

# --------------------------------------------------------------------------
# Step 1: LOO array (12 jobs, all run in parallel)
# --------------------------------------------------------------------------
echo ""
echo "=== Step 1: Submitting LOO array (12 jobs) ==="
LOO_JOB=$(sbatch --parsable scripts/cluster/run_loo.sh)
echo "  LOO job array ID: $LOO_JOB"

# --------------------------------------------------------------------------
# Step 2: Inference pass (runs after ALL 12 LOO jobs finish — not required
#         by inference itself, but keeps the pipeline sequential for monitoring)
# --------------------------------------------------------------------------
echo ""
echo "=== Step 2: Submitting inference pass (depends on LOO array) ==="
INFER_JOB=$(sbatch --parsable \
    --dependency=afterok:"$LOO_JOB" \
    scripts/cluster/run_inference.sh)
echo "  Inference job ID: $INFER_JOB"

# --------------------------------------------------------------------------
# Step 3: Analysis (depends on inference)
# --------------------------------------------------------------------------
echo ""
echo "=== Step 3: Submitting analysis (depends on inference) ==="
ANALYSIS_JOB=$(sbatch --parsable \
    --dependency=afterok:"$INFER_JOB" \
    scripts/cluster/run_analysis.sh)
echo "  Analysis job ID: $ANALYSIS_JOB"

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
echo ""
echo "============================================================"
echo "Jobs submitted:"
echo "  LOO array:  $LOO_JOB  (12 parallel tasks, ~24h each)"
echo "  Inference:  $INFER_JOB  (after LOO; ~12h)"
echo "  Analysis:   $ANALYSIS_JOB  (after inference; ~2h)"
echo ""
echo "Monitor:"
echo "  squeue -u \$USER"
echo "  tail -f cluster_logs/loo_0.out     # live output for task 0 (generic+node2vec)"
echo "  tail -f cluster_logs/inference.out  # inference progress"
echo ""
echo "When done, rsync results back to your local machine:"
echo "  rsync -avz cluster:$PROJECT_DIR/cluster_results/ ./cluster_results/"
echo "  rsync -avz cluster:$PROJECT_DIR/outputs_prediction/analysis/ ./outputs_prediction/analysis/"
echo "============================================================"
