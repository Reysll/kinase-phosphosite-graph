# src/io_utils.py
from __future__ import annotations

import os
import pandas as pd

"""
Read a table from disk with sensible defaults based on file extension.

Supported:
    - .csv  -> comma separated
    - .tsv  -> tab separated
    - .txt  -> tab separated (common in bio)
    - .xlsx/.xls -> Excel
    - .parquet -> Parquet (faster for large tables)
    - fallback -> try csv, then tsv
"""

def read_table(path: str) -> pd.DataFrame:
    ext = os.path.splitext(path)[1].lower()

    # Parquet (fast cached format)
    if ext == ".parquet":
        return pd.read_parquet(path)

    # Excel
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)

    # CSV (comma separated)
    if ext == ".csv":
        return pd.read_csv(path)

    # TSV / TXT
    if ext in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t")

    # Gzip compressed text (tsv.gz or txt.gz)
    if ext == ".gz":
        return pd.read_csv(path, sep="\t", compression="infer")

    # Fallback
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.read_csv(path, sep="\t")
