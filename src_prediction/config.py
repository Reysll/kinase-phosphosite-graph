from __future__ import annotations

from pathlib import Path


# Base project paths
PROJECT_ROOT = Path(".")
OUTPUTS_DIR = PROJECT_ROOT / "outputs"
PRED_OUTPUTS_DIR = PROJECT_ROOT / "outputs_prediction"

# Graph inputs from src/
GENERIC_NODES = OUTPUTS_DIR / "nodes.csv.gz"
GENERIC_EDGES = OUTPUTS_DIR / "edges.csv.gz"
LIVER_NODES = OUTPUTS_DIR / "Liver_network_nodes.csv.gz"
LIVER_EDGES = OUTPUTS_DIR / "Liver_network_edges.csv.gz"

# Prediction prep outputs
LIVER_SITE_NODES_OUT = PRED_OUTPUTS_DIR / "liver_site_nodes.csv.gz"
CANDIDATE_KINASES_OUT = PRED_OUTPUTS_DIR / "candidate_kinases.csv.gz"
LIVER_POSITIVE_EDGES_OUT = PRED_OUTPUTS_DIR / "liver_positive_kinase_site_edges.csv.gz"
PREP_SUMMARY_OUT = PRED_OUTPUTS_DIR / "prep_summary.txt"