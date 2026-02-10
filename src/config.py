# src/config.py
from __future__ import annotations
import os

# Project root = one folder above /src
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs")

# Input files (keep all paths consistent)
KINASE_SUBSTRATE_FILE = os.path.join(DATA_DIR, "Kinase_Substrate_Dataset")
PTMSIGDB_FILE = os.path.join(DATA_DIR, "PTMsigDB.txt")
PPASE_FILE = os.path.join(DATA_DIR, "PPase_protSubtrates_201903.xls")

MSIGDB_FILE = os.path.join(DATA_DIR, "MsigDB.txt")
PPI_FILE = os.path.join(DATA_DIR, "high_confidence_score(in).csv")
GENE_PROTEIN_MAP = os.path.join(DATA_DIR, "gene_protein(in).csv")
PTMCODE_FILE = os.path.join(DATA_DIR, "PTMcode2_associations_between_proteins.txt.gz")

# Safeguards (for later steps if needed)
MAX_SITE_CLIQUE = 50
MAX_PROTEIN_CLIQUE = 100
