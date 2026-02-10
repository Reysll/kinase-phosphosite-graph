# src/io_utils.py
from __future__ import annotations

import os
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    """
    Read a table from CSV/TSV or Excel (xls/xlsx).
    """
    ext = os.path.splitext(path)[1].lower()

    if ext in {".xls", ".xlsx"}:
        # Default: first sheet
        return pd.read_excel(path)

    # Try TSV then CSV for text files
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)
