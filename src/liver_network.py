from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, Set, Tuple

import numpy as np
import pandas as pd

from src.ids import protein_id, site_id
from src.io_utils import read_nodes_csv_gz, read_edges_csv_gz, write_nodes_csv_gz, write_edges_csv_gz

from src.builders.ptmsigdb import build_ptmsigdb_edges
from src.builders.msigdb import build_msigdb_edges
from src.builders.ppi import build_ppi_edges
from src.builders.ptmcode2 import build_ptmcode2_edges


@dataclass
class CorrEdge:
    source: str
    target: str
    relation: str
    weight: float
    n_common: int


PHOSPHO_SITE_RE = re.compile(r"([STY])\s*([0-9]+)")


def _merge_edges(base_edges: pd.DataFrame, new_edges: pd.DataFrame) -> pd.DataFrame:
    for col in ["source", "target", "relation"]:
        if col not in base_edges.columns:
            raise ValueError(f"base_edges missing required column: {col}")
        if col not in new_edges.columns:
            raise ValueError(f"new_edges missing required column: {col}")

    for col in base_edges.columns:
        if col not in new_edges.columns:
            new_edges[col] = np.nan
    for col in new_edges.columns:
        if col not in base_edges.columns:
            base_edges[col] = np.nan

    out = pd.concat([base_edges, new_edges], ignore_index=True)
    out = out.drop_duplicates(subset=["source", "target", "relation"])
    return out


def _ensure_node_metadata(nodes_df: pd.DataFrame) -> pd.DataFrame:
    out = nodes_df.copy()

    if "is_kinase" not in out.columns:
        out["is_kinase"] = np.nan
    if "is_phosphatase" not in out.columns:
        out["is_phosphatase"] = np.nan
    if "protein_role" not in out.columns:
        out["protein_role"] = np.nan

    protein_mask = out["node_type"] == "protein"
    site_mask = out["node_type"] == "site"

    out.loc[protein_mask & out["is_kinase"].isna(), "is_kinase"] = 0
    out.loc[protein_mask & out["is_phosphatase"].isna(), "is_phosphatase"] = 0
    out.loc[site_mask & out["is_kinase"].isna(), "is_kinase"] = 0
    out.loc[site_mask & out["is_phosphatase"].isna(), "is_phosphatase"] = 0

    out.loc[protein_mask & out["protein_role"].isna(), "protein_role"] = "protein"
    out.loc[site_mask & out["protein_role"].isna(), "protein_role"] = "site"

    out["is_kinase"] = out["is_kinase"].astype(int)
    out["is_phosphatase"] = out["is_phosphatase"].astype(int)

    return out


def _find_control_sample_columns(cols: List[str]) -> Tuple[List[str], List[str]]:
    control_cols: List[str] = []
    sample_cols: List[str] = []
    for c in cols:
        if isinstance(c, str) and c.strip().endswith("Control"):
            control_cols.append(c)
        elif isinstance(c, str) and c.strip().endswith("Sample"):
            sample_cols.append(c)
    return control_cols, sample_cols


def _identified_at_least(df: pd.DataFrame, cols: List[str], min_n: int) -> pd.Series:
    return df[cols].notna().sum(axis=1) >= min_n


def _extract_single_site_label(mod_str: str) -> Optional[str]:
    if not isinstance(mod_str, str):
        return None
    hits = PHOSPHO_SITE_RE.findall(mod_str)
    if len(hits) != 1:
        return None
    aa, pos = hits[0]
    return f"{aa}{pos}"


def _build_liver_additions(
    protein_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    min_identified: int = 6,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if "Gene Symbol" not in protein_df.columns:
        raise ValueError("ProteinExpression sheet missing required column: 'Gene Symbol'")

    prot_control, prot_sample = _find_control_sample_columns(list(protein_df.columns))
    prot_cols = prot_control + prot_sample
    if len(prot_cols) == 0:
        raise ValueError("ProteinExpression: could not detect any Control/Sample abundance columns.")

    keep_prot = _identified_at_least(protein_df, prot_cols, min_identified)
    prot_genes = (
        protein_df.loc[keep_prot, "Gene Symbol"]
        .astype(str)
        .str.strip()
        .replace({"": np.nan})
        .dropna()
        .unique()
        .tolist()
    )

    protein_nodes = [(protein_id(g), "protein") for g in prot_genes]

    if "Gene Symbol" not in phospho_df.columns:
        raise ValueError("Phosphorylation sheet missing required column: 'Gene Symbol'")
    if "Modifications in Master Proteins" not in phospho_df.columns:
        raise ValueError("Phosphorylation sheet missing required column: 'Modifications in Master Proteins'")

    ph_control, ph_sample = _find_control_sample_columns(list(phospho_df.columns))
    ph_cols = ph_control + ph_sample
    if len(ph_cols) == 0:
        raise ValueError("Phosphorylation: could not detect any Control/Sample abundance columns.")

    keep_ph = _identified_at_least(phospho_df, ph_cols, min_identified)
    ph_sub = phospho_df.loc[keep_ph, ["Gene Symbol", "Modifications in Master Proteins"] + ph_cols].copy()

    site_nodes: List[Tuple[str, str]] = []
    has_site_edges: List[Tuple[str, str, str]] = []

    for _, row in ph_sub.iterrows():
        gene = str(row["Gene Symbol"]).strip()
        if not gene:
            continue
        site_label = _extract_single_site_label(row["Modifications in Master Proteins"])
        if site_label is None:
            continue

        p = protein_id(gene)
        s = site_id(gene, site_label)
        site_nodes.append((s, "site"))
        has_site_edges.append((p, s, "has_site"))

    nodes_df = pd.DataFrame(protein_nodes + site_nodes, columns=["node_id", "node_type"]).drop_duplicates()
    edges_df = pd.DataFrame(has_site_edges, columns=["source", "target", "relation"]).drop_duplicates()

    return nodes_df, edges_df


def _prepare_site_fold_change_matrix(
    phospho_df: pd.DataFrame,
    min_identified: int = 6,
) -> Tuple[pd.DataFrame, List[str]]:
    if "Gene Symbol" not in phospho_df.columns:
        raise ValueError("Phosphorylation sheet missing required column: 'Gene Symbol'")
    if "Modifications in Master Proteins" not in phospho_df.columns:
        raise ValueError("Phosphorylation sheet missing required column: 'Modifications in Master Proteins'")

    ph_control, ph_sample = _find_control_sample_columns(list(phospho_df.columns))
    if len(ph_control) == 0 or len(ph_sample) == 0:
        raise ValueError("Phosphorylation: could not detect Control and Sample columns.")

    site_ids_all: List[str] = []
    for g, mod in zip(
        phospho_df["Gene Symbol"].astype(str).tolist(),
        phospho_df["Modifications in Master Proteins"].tolist(),
    ):
        gene = str(g).strip()
        site_label = _extract_single_site_label(mod)
        if gene and site_label:
            site_ids_all.append(site_id(gene, site_label))
        else:
            site_ids_all.append("")

    control_mat = phospho_df[ph_control].apply(pd.to_numeric, errors="coerce")
    sample_mat = phospho_df[ph_sample].apply(pd.to_numeric, errors="coerce")

    valid_site = pd.Series([sid != "" for sid in site_ids_all], index=phospho_df.index)
    keep = (
        valid_site
        & (control_mat.notna().sum(axis=1) >= 1)
        & (sample_mat.notna().sum(axis=1) >= min_identified)
    )

    control_mat = control_mat.loc[keep].copy()
    sample_mat = sample_mat.loc[keep].copy()
    kept_site_ids = [sid for sid, k in zip(site_ids_all, keep.tolist()) if k]

    mu = control_mat.mean(axis=1, skipna=True)
    mu = mu.where(mu > 0, np.nan)

    fold_change_mat = sample_mat.div(mu, axis=0)
    keep_fc = fold_change_mat.notna().sum(axis=1) >= min_identified

    fold_change_mat = fold_change_mat.loc[keep_fc].copy()
    final_site_ids = [sid for sid, k in zip(kept_site_ids, keep_fc.tolist()) if k]

    fold_change_mat.index = final_site_ids

    n_rows_before_collapse = len(fold_change_mat)
    n_unique_sites = fold_change_mat.index.nunique()

    if n_rows_before_collapse != n_unique_sites:
        fold_change_mat = fold_change_mat.groupby(level=0).mean()

    final_site_ids = fold_change_mat.index.tolist()

    print("=== Site fold-change matrix stats ===")
    print(f"Single parsed site rows kept: {n_rows_before_collapse:,}")
    print(f"Unique sites after duplicate collapse: {n_unique_sites:,}")
    print(f"Tumor sample columns used: {len(ph_sample):,}")

    return fold_change_mat, final_site_ids


def _prepare_protein_fold_change_matrix(
    protein_df: pd.DataFrame,
    min_identified: int = 6,
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Protein-level fold-change matrix:
    1. mean across healthy controls per protein
    2. tumor / mean_control per tumor sample
    3. collapse duplicate proteins by averaging rows
    """
    if "Gene Symbol" not in protein_df.columns:
        raise ValueError("ProteinExpression sheet missing required column: 'Gene Symbol'")

    prot_control, prot_sample = _find_control_sample_columns(list(protein_df.columns))
    if len(prot_control) == 0 or len(prot_sample) == 0:
        raise ValueError("ProteinExpression: could not detect Control and Sample columns.")

    protein_ids_all = []
    for g in protein_df["Gene Symbol"].astype(str).tolist():
        gene = str(g).strip()
        protein_ids_all.append(protein_id(gene) if gene else "")

    control_mat = protein_df[prot_control].apply(pd.to_numeric, errors="coerce")
    sample_mat = protein_df[prot_sample].apply(pd.to_numeric, errors="coerce")

    valid_protein = pd.Series([pid != "" for pid in protein_ids_all], index=protein_df.index)
    keep = (
        valid_protein
        & (control_mat.notna().sum(axis=1) >= 1)
        & (sample_mat.notna().sum(axis=1) >= min_identified)
    )

    control_mat = control_mat.loc[keep].copy()
    sample_mat = sample_mat.loc[keep].copy()
    kept_protein_ids = [pid for pid, k in zip(protein_ids_all, keep.tolist()) if k]

    mu = control_mat.mean(axis=1, skipna=True)
    mu = mu.where(mu > 0, np.nan)

    fold_change_mat = sample_mat.div(mu, axis=0)
    keep_fc = fold_change_mat.notna().sum(axis=1) >= min_identified

    fold_change_mat = fold_change_mat.loc[keep_fc].copy()
    final_protein_ids = [pid for pid, k in zip(kept_protein_ids, keep_fc.tolist()) if k]

    fold_change_mat.index = final_protein_ids

    n_rows_before_collapse = len(fold_change_mat)
    n_unique_proteins = fold_change_mat.index.nunique()

    if n_rows_before_collapse != n_unique_proteins:
        fold_change_mat = fold_change_mat.groupby(level=0).mean()

    final_protein_ids = fold_change_mat.index.tolist()

    print("=== Protein fold-change matrix stats ===")
    print(f"Protein rows kept: {n_rows_before_collapse:,}")
    print(f"Unique proteins after duplicate collapse: {n_unique_proteins:,}")
    print(f"Tumor sample columns used: {len(prot_sample):,}")

    return fold_change_mat, final_protein_ids


def _build_corr_edges_from_fc_matrix(
    fc_mat: pd.DataFrame,
    percentile: float,
    min_common: int,
    pos_relation: str,
    neg_relation: str,
    label: str,
) -> pd.DataFrame:
    if fc_mat.shape[0] < 2:
        print(f"Not enough {label} rows to compute correlations.")
        return pd.DataFrame(columns=["source", "target", "relation", "weight", "n_common"])

    corr = fc_mat.T.corr(method="pearson", min_periods=min_common)

    mask = np.triu(np.ones(corr.shape, dtype=bool), k=1)
    upper = corr.where(mask)
    corr_values = upper.stack().astype(float)

    positive_vals = corr_values[corr_values > 0]
    negative_vals = corr_values[corr_values < 0]

    pos_threshold = np.percentile(positive_vals, percentile) if len(positive_vals) > 0 else np.inf
    neg_threshold = np.percentile(negative_vals, 100 - percentile) if len(negative_vals) > 0 else -np.inf

    print(f"=== {label} fold-change correlation thresholds ===")
    print(f"Percentile: {percentile}")
    print(f"Positive threshold: {pos_threshold:.6f}" if np.isfinite(pos_threshold) else "Positive threshold: inf")
    print(f"Negative threshold: {neg_threshold:.6f}" if np.isfinite(neg_threshold) else "Negative threshold: -inf")
    print(f"Total pairwise correlations evaluated: {len(corr_values):,}")

    selected = upper.stack().reset_index()
    selected.columns = ["source", "target", "weight"]

    selected = selected.loc[
        (selected["weight"] >= pos_threshold) | (selected["weight"] <= neg_threshold)
    ].copy()

    if selected.empty:
        print(f"No {label} correlation edges passed threshold.")
        return pd.DataFrame(columns=["source", "target", "relation", "weight", "n_common"])

    selected["relation"] = np.where(
        selected["weight"] >= pos_threshold,
        pos_relation,
        neg_relation,
    )

    n_common_vals = []
    for row in selected.itertuples(index=False):
        common = int((fc_mat.loc[row.source].notna() & fc_mat.loc[row.target].notna()).sum())
        n_common_vals.append(common)

    selected["n_common"] = n_common_vals

    out = selected[["source", "target", "relation", "weight", "n_common"]].drop_duplicates(
        subset=["source", "target", "relation"]
    )

    print(f"{label} fold-change correlation edges kept: {len(out):,}")
    return out


def _build_site_fold_change_corr_edges(
    phospho_df: pd.DataFrame,
    min_identified: int = 6,
    min_common: int = 6,
    percentile: float = 80.0,
) -> pd.DataFrame:
    fc_mat, _ = _prepare_site_fold_change_matrix(
        phospho_df=phospho_df,
        min_identified=min_identified,
    )
    return _build_corr_edges_from_fc_matrix(
        fc_mat=fc_mat,
        percentile=percentile,
        min_common=min_common,
        pos_relation="site_corr_fc_pos",
        neg_relation="site_corr_fc_neg",
        label="Site",
    )


def _build_protein_fold_change_corr_edges(
    protein_df: pd.DataFrame,
    min_identified: int = 6,
    min_common: int = 6,
    percentile: float = 80.0,
) -> pd.DataFrame:
    fc_mat, _ = _prepare_protein_fold_change_matrix(
        protein_df=protein_df,
        min_identified=min_identified,
    )
    return _build_corr_edges_from_fc_matrix(
        fc_mat=fc_mat,
        percentile=percentile,
        min_common=min_common,
        pos_relation="protein_corr_fc_pos",
        neg_relation="protein_corr_fc_neg",
        label="Protein",
    )


def run(
    liver_xlsx: Path,
    protein_sheet: str,
    phospho_sheet: str,
    generic_nodes: Path,
    generic_edges: Path,
    out_nodes: Path,
    out_edges: Path,
    ptmsigdb_path: str,
    msigdb_path: str,
    ppi_path: str,
    gene_map_path: str,
    ptmcode_path: str,
    min_identified: int = 6,
    min_common: int = 6,
    site_corr_percentile: float = 80.0,
    protein_corr_percentile: float = 80.0,
    add_protein_fc_corr: bool = True,
    ppi_min_confidence: int = 0,
    ptmcode_chunksize: int = 250000,
) -> None:
    nodes_df = read_nodes_csv_gz(generic_nodes)
    edges_df = read_edges_csv_gz(generic_edges)

    if "weight" not in edges_df.columns:
        edges_df["weight"] = np.nan
    if "n_common" not in edges_df.columns:
        edges_df["n_common"] = np.nan

    protein_df = pd.read_excel(liver_xlsx, sheet_name=protein_sheet, header=1)
    phospho_df = pd.read_excel(liver_xlsx, sheet_name=phospho_sheet, header=1)

    new_nodes, new_has_site = _build_liver_additions(
        protein_df=protein_df,
        phospho_df=phospho_df,
        min_identified=min_identified,
    )

    nodes_before = nodes_df["node_id"].nunique()
    edges_before = edges_df.drop_duplicates(subset=["source", "target", "relation"]).shape[0]

    nodes_df = pd.concat([nodes_df, new_nodes], ignore_index=True).drop_duplicates(subset=["node_id", "node_type"])
    nodes_df = _ensure_node_metadata(nodes_df)

    edges_df = _merge_edges(edges_df, new_has_site.assign(weight=np.nan, n_common=np.nan))

    protein_ids: Set[str] = set(nodes_df.loc[nodes_df["node_type"] == "protein", "node_id"].astype(str))
    site_ids: Set[str] = set(nodes_df.loc[nodes_df["node_type"] == "site", "node_id"].astype(str))

    print("=== Reconnecting liver-added nodes via pathway/PPI/PTMcode datasets ===")

    ptm_edges = build_ptmsigdb_edges(ptmsigdb_path, allowed_site_ids=site_ids)
    edges_df = _merge_edges(edges_df, ptm_edges.assign(weight=np.nan, n_common=np.nan))

    msig_edges = build_msigdb_edges(msigdb_path, allowed_protein_ids=protein_ids)
    edges_df = _merge_edges(edges_df, msig_edges.assign(weight=np.nan, n_common=np.nan))

    ppi_edges = build_ppi_edges(
        ppi_path=ppi_path,
        gene_map_path=gene_map_path,
        existing_protein_ids=protein_ids,
        min_confidence_score=ppi_min_confidence,
    )
    if "confidence_score" in ppi_edges.columns:
        ppi_edges = ppi_edges.rename(columns={"confidence_score": "weight"})
    edges_df = _merge_edges(edges_df, ppi_edges)

    ptmcode_edges = build_ptmcode2_edges(
        ptmcode_path=ptmcode_path,
        gene_map_csv=gene_map_path,
        existing_site_ids=site_ids,
        chunksize=ptmcode_chunksize,
    )
    edges_df = _merge_edges(edges_df, ptmcode_edges.assign(weight=np.nan, n_common=np.nan))

    print("=== Adding liver site fold-change correlation edges ===")
    site_corr_edges = _build_site_fold_change_corr_edges(
        phospho_df=phospho_df,
        min_identified=min_identified,
        min_common=min_common,
        percentile=site_corr_percentile,
    )
    edges_df = _merge_edges(edges_df, site_corr_edges)

    if add_protein_fc_corr:
        print("=== Adding liver protein fold-change correlation edges ===")
        protein_corr_edges = _build_protein_fold_change_corr_edges(
            protein_df=protein_df,
            min_identified=min_identified,
            min_common=min_common,
            percentile=protein_corr_percentile,
        )
        edges_df = _merge_edges(edges_df, protein_corr_edges)

    write_nodes_csv_gz(nodes_df, out_nodes)
    write_edges_csv_gz(edges_df, out_edges)

    print("=== Liver network build complete ===")
    print(f"Generic nodes: {nodes_before:,}")
    print(f"Liver network nodes: {nodes_df['node_id'].nunique():,} (added {nodes_df['node_id'].nunique() - nodes_before:,})")
    print(f"Generic edges: {edges_before:,}")
    print(f"Liver network edges: {edges_df.drop_duplicates(subset=['source','target','relation']).shape[0]:,}")
    print(f"Wrote: {out_nodes}")
    print(f"Wrote: {out_edges}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--liver-xlsx", type=Path, required=True)
    ap.add_argument("--protein-sheet", type=str, default="ProteinExpression")
    ap.add_argument("--phospho-sheet", type=str, default="Phosphorylation")

    ap.add_argument("--generic-nodes", type=Path, default=Path("outputs/nodes.csv.gz"))
    ap.add_argument("--generic-edges", type=Path, default=Path("outputs/edges.csv.gz"))

    ap.add_argument("--out-nodes", type=Path, default=Path("outputs/Liver_network_nodes.csv.gz"))
    ap.add_argument("--out-edges", type=Path, default=Path("outputs/Liver_network_edges.csv.gz"))

    ap.add_argument("--ptmsigdb", type=str, default="data/PTMsigDB.txt")
    ap.add_argument("--msigdb", type=str, default="data/MsigDB.txt")
    ap.add_argument("--ppi", type=str, default="data/high_confidence_score.csv")
    ap.add_argument("--gene-map", type=str, default="data/gene_protein.csv")
    ap.add_argument("--ptmcode2", type=str, default="data/PTMcode2_associations_between_proteins.txt.gz")

    ap.add_argument("--min-identified", type=int, default=6)
    ap.add_argument("--min-common", type=int, default=6)
    ap.add_argument("--site-corr-percentile", type=float, default=80.0)
    ap.add_argument("--protein-corr-percentile", type=float, default=80.0)
    ap.add_argument("--add-protein-fc-corr", action="store_true")

    ap.add_argument("--ppi-min-confidence", type=int, default=0)
    ap.add_argument("--ptmcode-chunksize", type=int, default=250000)

    args = ap.parse_args()

    run(
        liver_xlsx=args.liver_xlsx,
        protein_sheet=args.protein_sheet,
        phospho_sheet=args.phospho_sheet,
        generic_nodes=args.generic_nodes,
        generic_edges=args.generic_edges,
        out_nodes=args.out_nodes,
        out_edges=args.out_edges,
        ptmsigdb_path=args.ptmsigdb,
        msigdb_path=args.msigdb,
        ppi_path=args.ppi,
        gene_map_path=args.gene_map,
        ptmcode_path=args.ptmcode2,
        min_identified=args.min_identified,
        min_common=args.min_common,
        site_corr_percentile=args.site_corr_percentile,
        protein_corr_percentile=args.protein_corr_percentile,
        add_protein_fc_corr=args.add_protein_fc_corr,
        ppi_min_confidence=args.ppi_min_confidence,
        ptmcode_chunksize=args.ptmcode_chunksize,
    )


if __name__ == "__main__":
    main()