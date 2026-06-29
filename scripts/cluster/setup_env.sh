#!/usr/bin/env bash
# =============================================================================
# Environment setup for UTRGV CRADLE HPC cluster (login.cradle.utrgv.edu).
# Run this ONCE on the login node before submitting any jobs.
#
# Prerequisites (one-time, run manually before this script):
#   curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
#   bash Miniconda3-latest-Linux-x86_64.sh -b -p ~/miniconda3
#   ~/miniconda3/bin/conda create -n kinase python=3.12 -y
#
# Then run this script:
#   cd ~/Kinase_Substrate_HGN
#   bash scripts/cluster/setup_env.sh
#
# The SLURM job scripts activate cluster_venv/ automatically.
# =============================================================================

set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
VENV_DIR="$PROJECT_DIR/cluster_venv"
CONDA_PYTHON="$HOME/miniconda3/envs/kinase/bin/python"

echo "=== Project root: $PROJECT_DIR"
echo "=== Venv target:  $VENV_DIR"

# --------------------------------------------------------------------------
# Verify Python 3.12 from the conda environment
# System python3 is 3.6.8 (too old); we use conda's Python 3.12 instead.
# --------------------------------------------------------------------------
if [ ! -f "$CONDA_PYTHON" ]; then
    echo "ERROR: $CONDA_PYTHON not found."
    echo "Run the prerequisite commands in the script header first."
    exit 1
fi

echo "=== Using Python: $CONDA_PYTHON ==="
"$CONDA_PYTHON" --version

# --------------------------------------------------------------------------
# Create virtual environment using Python 3.12
# --------------------------------------------------------------------------
echo "=== Creating virtual environment at $VENV_DIR ==="
"$CONDA_PYTHON" -m venv "$VENV_DIR"
source "$VENV_DIR/bin/activate"

python --version
pip install --upgrade pip --quiet

# --------------------------------------------------------------------------
# Install dependencies
# --------------------------------------------------------------------------
echo "=== Installing dependencies from requirements.txt ==="
pip install -r "$PROJECT_DIR/requirements.txt"

echo ""
echo "=== Setup complete. ==="
echo "Venv: $VENV_DIR"
echo "Test: source $VENV_DIR/bin/activate && python -c 'import networkx, gensim, sklearn; print(\"OK\")'"
