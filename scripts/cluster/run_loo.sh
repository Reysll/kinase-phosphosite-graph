#!/usr/bin/env bash
# =============================================================================
# SLURM job array: LOO evaluation, all 12 (network x embedding) combos.
#
# Each array task runs ONE (network, embedding) combo independently.
# All 12 run in parallel (subject to the cluster's scheduling; kimq has 2 nodes
# so at most 2 tasks run at once — the rest queue).
#
# Submit: sbatch scripts/cluster/run_loo.sh
# Monitor: squeue -u $USER
# Logs: cluster_logs/loo_<task_id>.{out,err}
# =============================================================================

# Cluster specs (UTRGV CRADLE login.cradle.utrgv.edu, as of 2026-06-28):
#   Partition: kimq  |  2 nodes (sxm003, sxm004)
#   CPUs/node: 56 (2 reserved for OS → 54 usable; requesting 48 to be safe)
#   RAM/node:  1 TB
#   Wall limit: 12 h

#SBATCH --job-name=kinase_loo
#SBATCH --array=0-11                  # 12 combos: 4 networks x 3 embeddings
#SBATCH --partition=kimq
#SBATCH --cpus-per-task=48
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --output=cluster_logs/loo_%a.out
#SBATCH --error=cluster_logs/loo_%a.err

set -euo pipefail

# --------------------------------------------------------------------------
# CPU allocation: 36 outer trial workers + 12 node2vec Word2Vec threads = 48
# --------------------------------------------------------------------------
N_JOBS=36
WORKERS=12
DIMENSIONS=64
TRIALS=all          # "all" = all sites (single + multi-kinase); "multi" = original 6,909
OUTPUT_DIR=cluster_results

# --------------------------------------------------------------------------
# Task ID → (network, embedding) mapping
# 0..2  → generic    + node2vec / spectral / concat
# 3..5  → control    + node2vec / spectral / concat
# 6..8  → cancer     + node2vec / spectral / concat
# 9..11 → liver_fc   + node2vec / spectral / concat
# --------------------------------------------------------------------------
NETWORKS=(generic control cancer liver_fc)
EMBEDDINGS=(node2vec spectral concat)

TASK_ID=${SLURM_ARRAY_TASK_ID:-0}
NETWORK=${NETWORKS[$((TASK_ID / 3))]}
EMBEDDING=${EMBEDDINGS[$((TASK_ID % 3))]}

echo "=== SLURM task $TASK_ID: $NETWORK x $EMBEDDING ==="
echo "n_jobs=$N_JOBS  workers=$WORKERS  dimensions=$DIMENSIONS  trials=$TRIALS"
date

# --------------------------------------------------------------------------
# Environment
# --------------------------------------------------------------------------
PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VENV_DIR="$PROJECT_DIR/cluster_venv"

source "$VENV_DIR/bin/activate"
cd "$PROJECT_DIR"
export PYTHONPATH="$PROJECT_DIR:${PYTHONPATH:-}"

mkdir -p cluster_logs "$OUTPUT_DIR"

# --------------------------------------------------------------------------
# Run
# --------------------------------------------------------------------------
python scripts/cluster/cluster_loo_runner.py \
    --network    "$NETWORK"    \
    --embedding  "$EMBEDDING"  \
    --dimensions "$DIMENSIONS" \
    --n-jobs     "$N_JOBS"     \
    --workers    "$WORKERS"    \
    --trials     "$TRIALS"     \
    --output-dir "$OUTPUT_DIR"

echo "=== Task $TASK_ID done ==="
date
