from __future__ import annotations

from pathlib import Path
import os
import pandas as pd


def read_table(path: str | Path) -> pd.DataFrame:
    """
    Generic table reader based on file extension.
    Supports parquet, excel, csv, tsv/txt, and gz.
    """
    path = str(path)
    ext = os.path.splitext(path)[1].lower()

    # Parquet
    if ext == ".parquet":
        return pd.read_parquet(path)

    # Excel
    if ext in [".xlsx", ".xls"]:
        return pd.read_excel(path)

    # CSV (compressed or not)
    if ext == ".csv":
        return pd.read_csv(path, compression="infer")

    # TSV / TXT (compressed or not)
    if ext in [".tsv", ".txt"]:
        return pd.read_csv(path, sep="\t", compression="infer")

    # .gz fallback
    if ext == ".gz":
        if path.endswith(".csv.gz") or ".csv" in path:
            return pd.read_csv(path, compression="infer")
        return pd.read_csv(path, sep="\t", compression="infer")

    raise ValueError(f"Unsupported file format: {path}")


def read_nodes_csv_gz(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, compression="gzip")


def read_edges_csv_gz(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, compression="gzip")


def write_nodes_csv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, compression="gzip")


def write_edges_csv_gz(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False, compression="gzip")
