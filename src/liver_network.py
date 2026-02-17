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


PHOSPHO_SITE_RE = re.compile(r"([STY])\s*([0-9]+)")  # catches S311, T539, Y476


def _merge_edges(base_edges: pd.DataFrame, new_edges: pd.DataFrame) -> pd.DataFrame:
    """
    Your project rule: deduplicate on (source, target, relation).
    base_edges may optionally have weight/n_common; we preserve extra columns.
    """
    # Ensure core columns exist
    for col in ["source", "target", "relation"]:
        if col not in base_edges.columns:
            raise ValueError(f"base_edges missing required column: {col}")
        if col not in new_edges.columns:
            raise ValueError(f"new_edges missing required column: {col}")

    # Align columns
    for col in base_edges.columns:
        if col not in new_edges.columns:
            new_edges[col] = np.nan
    for col in new_edges.columns:
        if col not in base_edges.columns:
            base_edges[col] = np.nan

    out = pd.concat([base_edges, new_edges], ignore_index=True)
    out = out.drop_duplicates(subset=["source", "target", "relation"])
    return out


def _find_control_sample_columns(cols: List[str]) -> Tuple[List[str], List[str]]:
    """
    Uses suffix naming in your Excel:
    - endswith('Control')
    - endswith('Sample')
    """
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
    """
    Example:
      'P02768 1xPhospho [S311(100)]' -> 'S311'
      '... [Y476(100); S478(100)]' -> None (multi-site, skip)
    """
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
    """
    Builds:
      - new nodes from liver sheets (proteins + single-site phosphosites)
      - structural has_site edges for those phosphosites
    Filtering rule: keep rows identified in >= min_identified across ALL samples (control+sample).
    """
    # ProteinExpression
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

    # Phosphorylation
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
            continue  # skip multi-site or unparseable

        p = protein_id(gene)
        s = site_id(gene, site_label)
        site_nodes.append((s, "site"))
        has_site_edges.append((p, s, "has_site"))

    nodes_df = pd.DataFrame(protein_nodes + site_nodes, columns=["node_id", "node_type"]).drop_duplicates()
    edges_df = pd.DataFrame(has_site_edges, columns=["source", "target", "relation"]).drop_duplicates()

    return nodes_df, edges_df


def _prep_matrix_for_corr(df: pd.DataFrame, id_list: List[str], cols: List[str], min_identified: int) -> Tuple[pd.DataFrame, List[str]]:
    """
    Returns numeric matrix with rows filtered by >=min_identified in the given cols.
    """
    mat = df[cols].apply(pd.to_numeric, errors="coerce")
    keep = mat.notna().sum(axis=1) >= min_identified
    mat = mat.loc[keep].copy()
    ids = [id_list[i] for i, k in enumerate(keep.tolist()) if k]
    return mat, ids


def _corr_edges_from_matrix(
    mat: pd.DataFrame,
    ids: List[str],
    relation: str,
    min_common: int = 6,
    abs_threshold: float = 0.0,
    top_k: int = 20,
) -> List[CorrEdge]:
    """
    Top-K correlation edges per node.
    For each node i, keep the K strongest correlations by absolute value
    (optionally also require abs(r) >= abs_threshold).
    Adds directed edges (i->j) for kept neighbors.
    """
    if mat.shape[0] < 2:
        return []

    corr = mat.T.corr(method="pearson", min_periods=min_common)
    vals = corr.values
    n = vals.shape[0]

    edges: List[CorrEdge] = []

    for i in range(n):
        r_row = vals[i, :].copy()
        r_row[i] = np.nan

        # Candidates that are not NaN and pass threshold
        candidates = np.where(~np.isnan(r_row) & (np.abs(r_row) >= abs_threshold))[0]
        if candidates.size == 0:
            continue

        # Sort candidates by abs correlation descending and take top_k
        order = candidates[np.argsort(-np.abs(r_row[candidates]))]
        if top_k is not None and top_k > 0:
            order = order[:top_k]

        xi = mat.iloc[i, :]
        u = ids[i]

        for j in order:
            xj = mat.iloc[j, :]
            common = int((xi.notna() & xj.notna()).sum())
            if common < min_common:
                continue
            v = ids[j]
            edges.append(CorrEdge(u, v, relation, float(r_row[j]), common))

    return edges



def _build_liver_correlation_edges(
    protein_df: pd.DataFrame,
    phospho_df: pd.DataFrame,
    min_identified: int = 6,
    min_common: int = 6,
    abs_threshold: float = 0.7,
) -> pd.DataFrame:
    """
    Builds correlation edges per Dr. Ayati rules:
    - control correlations
    - sample correlations
    - only compute if overlap >= 6
    - only keep if abs(r) >= abs_threshold (to keep graph bounded)
    """
    edges: List[CorrEdge] = []

    # Protein correlations
    prot_control, prot_sample = _find_control_sample_columns(list(protein_df.columns))
    if "Gene Symbol" in protein_df.columns and len(prot_control) > 0 and len(prot_sample) > 0:
        prot_ids_all = [protein_id(str(g).strip()) if str(g).strip() else "" for g in protein_df["Gene Symbol"].astype(str).tolist()]

        mat_c, ids_c = _prep_matrix_for_corr(protein_df, prot_ids_all, prot_control, min_identified)
        mat_s, ids_s = _prep_matrix_for_corr(protein_df, prot_ids_all, prot_sample, min_identified)

        edges += _corr_edges_from_matrix(mat_c, ids_c, "protein_corr_control", min_common, abs_threshold)
        edges += _corr_edges_from_matrix(mat_s, ids_s, "protein_corr_sample", min_common, abs_threshold)

    # Site correlations
    ph_control, ph_sample = _find_control_sample_columns(list(phospho_df.columns))
    if "Gene Symbol" in phospho_df.columns and "Modifications in Master Proteins" in phospho_df.columns and len(ph_control) > 0 and len(ph_sample) > 0:
        site_ids_all: List[str] = []
        for g, mod in zip(phospho_df["Gene Symbol"].astype(str).tolist(), phospho_df["Modifications in Master Proteins"].tolist()):
            gene = str(g).strip()
            lab = _extract_single_site_label(mod)
            if gene and lab:
                site_ids_all.append(site_id(gene, lab))
            else:
                site_ids_all.append("")

        # Filter invalid IDs by turning those rows into all-NaN so they get dropped by min_identified
        ph_tmp = phospho_df.copy()
        invalid = [i for i, sid in enumerate(site_ids_all) if sid == ""]
        if invalid:
            ph_tmp.loc[ph_tmp.index[invalid], ph_control + ph_sample] = np.nan

        mat_c, ids_c = _prep_matrix_for_corr(ph_tmp, site_ids_all, ph_control, min_identified)
        mat_s, ids_s = _prep_matrix_for_corr(ph_tmp, site_ids_all, ph_sample, min_identified)

        edges += _corr_edges_from_matrix(mat_c, ids_c, "site_corr_control", min_common, abs_threshold)
        edges += _corr_edges_from_matrix(mat_s, ids_s, "site_corr_sample", min_common, abs_threshold)

    out = pd.DataFrame(
        [(e.source, e.target, e.relation, e.weight, e.n_common) for e in edges],
        columns=["source", "target", "relation", "weight", "n_common"],
    ).drop_duplicates(subset=["source", "target", "relation"])

    return out


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
    abs_threshold: float = 0.7,
    ppi_min_confidence: int = 0,
    ptmcode_chunksize: int = 250000,
) -> None:
    # Load generic outputs
    nodes_df = read_nodes_csv_gz(generic_nodes)
    edges_df = read_edges_csv_gz(generic_edges)

    # Ensure optional columns exist for later merges
    if "weight" not in edges_df.columns:
        edges_df["weight"] = np.nan
    if "n_common" not in edges_df.columns:
        edges_df["n_common"] = np.nan

    # Read liver excel (header row 2 => header=1)
    protein_df = pd.read_excel(liver_xlsx, sheet_name=protein_sheet, header=1)
    phospho_df = pd.read_excel(liver_xlsx, sheet_name=phospho_sheet, header=1)

    # Build liver additions
    new_nodes, new_has_site = _build_liver_additions(
        protein_df=protein_df,
        phospho_df=phospho_df,
        min_identified=min_identified,
    )

    # Merge nodes + has_site
    nodes_before = nodes_df["node_id"].nunique()
    edges_before = edges_df.drop_duplicates(subset=["source", "target", "relation"]).shape[0]

    nodes_df = pd.concat([nodes_df, new_nodes], ignore_index=True).drop_duplicates(subset=["node_id", "node_type"])
    edges_df = _merge_edges(edges_df, new_has_site.assign(weight=np.nan, n_common=np.nan))

    # Recompute allowed sets
    protein_ids: Set[str] = set(nodes_df.loc[nodes_df["node_type"] == "protein", "node_id"].astype(str))
    site_ids: Set[str] = set(nodes_df.loc[nodes_df["node_type"] == "site", "node_id"].astype(str))

    # Re-run secondary datasets (builders return DataFrames with source/target/relation)
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
    # If PPI builder includes a confidence_score column, keep it as weight
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

    # Correlation edges
    print("=== Adding liver correlation edges ===")
    corr_edges = _build_liver_correlation_edges(
        protein_df=protein_df,
        phospho_df=phospho_df,
        min_identified=min_identified,
        min_common=min_common,
        abs_threshold=abs_threshold,
    )
    edges_df = _merge_edges(edges_df, corr_edges)

    # Write outputs
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
    ap.add_argument("--abs-threshold", type=float, default=0.9)

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
        abs_threshold=args.abs_threshold,
        ppi_min_confidence=args.ppi_min_confidence,
        ptmcode_chunksize=args.ptmcode_chunksize,
    )


if __name__ == "__main__":
    main()
