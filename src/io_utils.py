# src/io_utils.py
from __future__ import annotations

import os
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    """
    Read a table from disk with sensible defaults based on file extension.

    Supported:
      - .csv  -> comma separated
      - .tsv  -> tab separated
      - .txt  -> tab separated (common in bio)
      - .xlsx/.xls -> Excel
      - fallback -> try csv, then tsv
    """
    ext = os.path.splitext(path)[1].lower()

    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)

    if ext == ".csv":
        # CSV should be comma by default
        return pd.read_csv(path)

    if ext in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t")

    # fallback: try csv then tsv
    try:
        return pd.read_csv(path)
    except Exception:
        return pd.read_csv(path, sep="\t")
