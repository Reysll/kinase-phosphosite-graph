# config.py
# src/config.py
from __future__ import annotations
import os

# Project root = one folder above /src
PROJECT_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))

DATA_DIR = os.path.join(PROJECT_ROOT, "data")
OUTPUT_DIR = os.path.join(PROJECT_ROOT, "outputs")

# Input files
KINASE_SUBSTRATE_FILE = os.path.join(DATA_DIR, "Kinase_Substrate_Dataset")

# If your file has an extension, change it like:
# KINASE_SUBSTRATE_FILE = os.path.join(DATA_DIR, "Kinase_Substrate_Dataset.tsv")

DATA_DIR = "data"

KINASE_SUBSTRATE_FILE = f"{DATA_DIR}/Kinase_Substrate_Dataset"
PHOSPHATASE_FILE = f"{DATA_DIR}/PPase_protSubtrates_201903.xls"
PTMSIG_FILE = f"{DATA_DIR}/PTMsigDB.txt"
MSIG_FILE = f"{DATA_DIR}/MsigDB.txt"
PPI_FILE = f"{DATA_DIR}/high_confidence_score(in).csv"
GENE_PROTEIN_MAP = f"{DATA_DIR}/gene_protein(in).csv"
PTMCODE_FILE = f"{DATA_DIR}/PTMcode2_associations_between_proteins.txt.gz"

# Safeguards against edge explosion
MAX_SITE_CLIQUE = 50
MAX_PROTEIN_CLIQUE = 100

OUTPUT_DIR = "outputs"
