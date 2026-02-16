# src/preprocess.py
from __future__ import annotations

from pathlib import Path
import hashlib
import pandas as pd

from src.config import (
    DATA_DIR,
    KINASE_SUBSTRATE_FILE,
    PROCESSED_DIR,
)


# Ensure processed dir exists
Path(PROCESSED_DIR).mkdir(parents=True, exist_ok=True)


def _hash_file(path: Path, block_size: int = 1_048_576) -> str:
    h = hashlib.md5()
    with path.open("rb") as f:
        while True:
            b = f.read(block_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()[:10]


def _sniff_sep(path: Path) -> str:
    with path.open("rb") as f:
        sample = f.read(4096)
    text = sample.decode("utf-8", errors="ignore")
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.lstrip().startswith("#"):
            continue
        comma = line.count(",")
        tab = line.count("\t")
        return "\t" if tab > comma else ","
    return ","


def preprocess_kinase_substrate() -> str:
    """
    Build or reuse the processed human-only Kinase_Substrate_Dataset TSV.
    Returns path to processed file.
    """
    raw_path = Path(KINASE_SUBSTRATE_FILE)

    out_path = Path(PROCESSED_DIR) / "Kinase_Substrate_Dataset_human.tsv"

    if out_path.exists():
        print(f"[cache] Using existing processed file: {out_path}")
        return str(out_path)

    if not raw_path.exists():
        raise FileNotFoundError(f"Raw Kinase_Substrate file not found: {raw_path}")

    print(f"[build] Reading raw TSV: {raw_path}")

    df = pd.read_csv(raw_path, sep="\t")

    # Keep only human kinase + human substrate
    # These column names match what we used earlier:
    #   col D = kinase organism
    #   col I = substrate organism
    # If your headers differ, adjust here.
    # Safe fallback: use positional columns if needed.
    if "KIN_ORGANISM" in df.columns and "SUB_ORGANISM" in df.columns:
        kin_org_col = "KIN_ORGANISM"
        sub_org_col = "SUB_ORGANISM"
    else:
        # positional: D is index 3, I is index 8
        kin_org_col = df.columns[3]
        sub_org_col = df.columns[8]

    df_human = df[(df[kin_org_col] == "human") & (df[sub_org_col] == "human")].copy()

    df_human.to_csv(out_path, sep="\t", index=False)
    print(f"[done] Filtered rows: {len(df_human):,}")
    print(f"[write] Wrote processed TSV: {out_path}")

    return str(out_path)


def preprocess_ppi(ppi_path: str, min_confidence_score: int = 0) -> str:
    """
    Preprocess PPI into a clean 3-column gz CSV:
      protein1, protein2, confidence_score

    Output: data/processed/ppi_min{min_confidence_score}.csv.gz
    """
    ppi_path = Path(ppi_path)
    out_path = Path(PROCESSED_DIR) / f"ppi_min{min_confidence_score}.csv.gz"

    if out_path.exists():
        print(f"[cache] Using existing processed PPI: {out_path}")
        return str(out_path)

    if not ppi_path.exists():
        raise FileNotFoundError(f"PPI file not found: {ppi_path}")

    sep = _sniff_sep(ppi_path)
    print(f"[build] Reading PPI: {ppi_path} (sep={repr(sep)})")

    df = pd.read_csv(ppi_path, sep=sep, compression="infer", low_memory=False)

    expected = {"protein1", "protein2", "confidence_score"}
    if not expected.issubset(set(df.columns)):
        if df.shape[1] < 3:
            raise ValueError(f"PPI parsed with <3 columns. Columns={list(df.columns)[:10]}")
        df = df.iloc[:, :3].copy()
        df.columns = ["protein1", "protein2", "confidence_score"]
    else:
        df = df[["protein1", "protein2", "confidence_score"]].copy()

    df["confidence_score"] = pd.to_numeric(df["confidence_score"], errors="coerce").fillna(0).astype(int)
    df = df[df["confidence_score"] >= min_confidence_score].copy()

    df.to_csv(out_path, index=False, compression="gzip")
    print(f"[done] Wrote processed PPI: {out_path} rows={len(df):,}")
    return str(out_path)


def preprocess_ptmcode2(ptmcode_path: str) -> str:
    """
    Preprocess PTMcode2 into a smaller gzipped TSV while preserving schema expected
    by src/builders/ptmcode2.py.

    Output: data/processed/ptmcode2_filtered_<hash>.tsv.gz
    """
    ptmcode_path = Path(ptmcode_path)
    if not ptmcode_path.exists():
        raise FileNotFoundError(f"PTMcode2 file not found: {ptmcode_path}")

    sig = _hash_file(ptmcode_path)
    out_path = Path(PROCESSED_DIR) / f"ptmcode2_filtered_{sig}.tsv.gz"

    if out_path.exists():
        print(f"[cache] Using existing processed PTMcode2: {out_path}")
        return str(out_path)

    print(f"[build] Reading PTMcode2 (one-time slow step): {ptmcode_path}")

    cols = [
        "Protein1", "Protein2", "Species",
        "PTM1", "Residue1", "rRCS1", "Propagated1",
        "PTM2", "Residue2", "rRCS2", "Propagated2",
        "Coevolution_evidence", "Manual_evidence", "Structure_distance_evidence",
    ]

    first_write = True
    total_kept = 0

    for chunk in pd.read_csv(
        ptmcode_path,
        sep="\t",
        compression="infer",
        comment="#",
        header=None,
        names=cols,
        chunksize=500_000,
        low_memory=False,
    ):
        species = chunk["Species"].astype(str)
        ptm1 = chunk["PTM1"].astype(str).str.lower()
        ptm2 = chunk["PTM2"].astype(str).str.lower()
        coev = pd.to_numeric(chunk["Coevolution_evidence"], errors="coerce").fillna(0).astype(int)

        mask = (
            (species == "Homo sapiens")
            & (ptm1 == "phosphorylation")
            & (ptm2 == "phosphorylation")
            & (coev >= 1)
        )

        keep = chunk.loc[mask, cols].copy()
        if keep.empty:
            continue

        total_kept += len(keep)

        keep.to_csv(
            out_path,
            sep="\t",
            index=False,
            header=False,
            mode="wt" if first_write else "at",
            compression="gzip",
        )
        first_write = False

    print(f"[done] Wrote processed PTMcode2: {out_path} rows={total_kept:,}")
    return str(out_path)
