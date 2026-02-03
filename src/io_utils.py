# src/io_utils.py
import pandas as pd


def read_table(path: str) -> pd.DataFrame:
    """
    Try reading as TSV first, then CSV.
    Many biology files are tab-delimited even when unnamed.
    """
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        return pd.read_csv(path)
