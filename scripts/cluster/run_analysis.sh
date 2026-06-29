#!/usr/bin/env bash
# =============================================================================
# SLURM job: run all analysis scripts on cluster results.
#
# NOTE: Analysis scripts are fast (~minutes). The easier path is to rsync
# cluster_results/ back to your local machine and run analysis there.
# This script is provided for fully automated cluster-side analysis.
#
# Submit after inference:
#   sbatch --dependency=afterok:<INFERENCE_JOB_ID> scripts/cluster/run_analysis.sh
# (submit_all.sh handles this automatically)
# =============================================================================

#SBATCH --job-name=kinase_analysis
#SBATCH --partition=kimq
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --time=02:00:00
#SBATCH --output=cluster_logs/analysis.out
#SBATCH --error=cluster_logs/analysis.err

set -euo pipefail

echo "=== Analysis pass ==="
date

PROJECT_DIR="${SLURM_SUBMIT_DIR}"
VENV_DIR="$PROJECT_DIR/cluster_venv"

source "$VENV_DIR/bin/activate"
cd "$PROJECT_DIR"
export PYTHONPATH="$PROJECT_DIR:${PYTHONPATH:-}"

mkdir -p cluster_logs cluster_results/analysis

# -----------------------------------------------------------------------
# Symlink each cluster LOO experiment dir into outputs_prediction/ so
# analyze_results.py finds them without code changes.
# -----------------------------------------------------------------------
echo "=== Symlinking cluster LOO results into outputs_prediction/ ==="
for exp_dir in cluster_results/*/; do
    exp_name=$(basename "$exp_dir")
    [[ "$exp_name" == "inference_all" ]] && continue
    target="$PROJECT_DIR/outputs_prediction/$exp_name"
    if [ ! -e "$target" ]; then
        ln -s "$PROJECT_DIR/$exp_dir" "$target"
        echo "  Linked: $exp_name"
    fi
done

if [ -d cluster_results/inference_all ] && [ ! -e outputs_prediction/inference_all_cluster ]; then
    ln -s "$PROJECT_DIR/cluster_results/inference_all" \
          "$PROJECT_DIR/outputs_prediction/inference_all_cluster"
    echo "Linked inference_all"
fi

echo ""
echo "=== Running analyze_results ==="
python -m src_prediction.analysis.analyze_results

echo ""
echo "=== Running build_rank_heatmap ==="
python -m src_prediction.analysis.build_rank_heatmap

echo ""
echo "=== Running build_rank_ks_test ==="
python -m src_prediction.analysis.build_rank_ks_test

echo ""
echo "=== Running spotlight + volcano + ksea scatter ==="
python -m src_prediction.analysis.build_spotlight_plot
python -m src_prediction.analysis.build_ks_volcano
python -m src_prediction.analysis.build_ksea_scatter

echo ""
echo "=== Analysis done ==="
date
