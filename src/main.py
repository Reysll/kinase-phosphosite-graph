# src/main.py
from __future__ import annotations

import os
import pandas as pd

from src.builders.kinase_substrate import build_kinase_substrate_graph
from src.builders.ppase_substrate import build_ppase_substrate_graph
from src.builders.ptmsigdb import build_ptmsigdb_edges
from src.builders.msigdb import build_msigdb_edges
from src.builders.ppi import build_ppi_edges
from src.builders.ptmcode2 import build_ptmcode2_edges

from src.preprocess import preprocess_kinase_substrate, preprocess_ppi, preprocess_ptmcode2


from src.sanity_checks import check_site_collisions

from src.config import (
    PPASE_FILE,
    PTMSIGDB_FILE,
    MSIGDB_FILE,
    PPI_FILE,
    GENE_PROTEIN_MAP,
    PTMCODE_FILE,
    OUTPUT_DIR,
)


def _merge_nodes(nodes_a: pd.DataFrame, nodes_b: pd.DataFrame) -> pd.DataFrame:
    return (
        pd.concat([nodes_a, nodes_b], ignore_index=True)
        .drop_duplicates(subset=["node_id", "node_type"])
        .reset_index(drop=True)
    )


def _merge_edges(edges_a: pd.DataFrame, edges_b: pd.DataFrame) -> pd.DataFrame:
    """
    Merge edges and deduplicate.

    Important rule for edges that carry a confidence_score (PPI):
      - If multiple rows share the same (source, target, relation),
        keep the row with the MAX confidence_score.
    """
    merged = pd.concat([edges_a, edges_b], ignore_index=True)

    for col in ["source", "target", "relation"]:
        if col not in merged.columns:
            raise ValueError(f"Edges missing required column: {col}")

    triple = ["source", "target", "relation"]

    if "confidence_score" in merged.columns:
        tmp = merged.copy()
        tmp["_conf_tmp"] = tmp["confidence_score"].fillna(-1)
        idx = tmp.groupby(triple)["_conf_tmp"].idxmax()
        merged = merged.loc[idx].drop(columns=["_conf_tmp"], errors="ignore").reset_index(drop=True)
    else:
        merged = merged.drop_duplicates(subset=triple).reset_index(drop=True)

    return merged


def main() -> None:
    # 1) Kinase-substrate (cached filtered file)
    processed_kinase_file = preprocess_kinase_substrate()

    print("=== Building base graph from Kinase_Substrate_Dataset (processed human-only) ===")
    k_nodes, k_edges = build_kinase_substrate_graph(processed_kinase_file)
    print(f"Kinase_Substrate produced: nodes={len(k_nodes):,} edges={len(k_edges):,}\n")

    # 2) PPase
    print("=== Building base graph from PPase_protSubtrates ===")
    p_nodes, p_edges = build_ppase_substrate_graph(PPASE_FILE)
    print(f"PPase produced: nodes={len(p_nodes):,} edges={len(p_edges):,}\n")

    # Merge base graphs
    print("=== Merging base graphs (Kinase + PPase) ===")
    nodes_df = _merge_nodes(k_nodes, p_nodes)
    edges_df = _merge_edges(k_edges, p_edges)
    print(f"Base nodes total (unique node_id+type): {len(nodes_df):,}")
    print(f"Base edges total (unique triples):       {len(edges_df):,}\n")

    # Sanity checks
    check_site_collisions(nodes_df)
    print()

    # 3) PTMsigDB site-site edges
    print("=== Adding PTMsigDB site-site pathway edges (restricted to existing sites) ===")
    site_ids = set(nodes_df.loc[nodes_df["node_type"] == "site", "node_id"])
    ptm_edges = build_ptmsigdb_edges(PTMSIGDB_FILE, site_ids)
    before = len(edges_df)
    edges_df = _merge_edges(edges_df, ptm_edges)
    print(f"PTMsigDB produced edges:                {len(ptm_edges):,}")
    print(f"New PTMsigDB edges added after merge:   {len(edges_df) - before:,}")
    print(f"Total edges now (unique triples):       {len(edges_df):,}\n")

    # 4) MsigDB protein-protein edges
    print("=== Adding MsigDB protein-protein pathway edges (restricted to existing proteins) ===")
    protein_ids = set(nodes_df.loc[nodes_df["node_type"] == "protein", "node_id"])
    msig_edges = build_msigdb_edges(MSIGDB_FILE, protein_ids)
    before = len(edges_df)
    edges_df = _merge_edges(edges_df, msig_edges)
    print(f"MsigDB produced edges:                  {len(msig_edges):,}")
    print(f"New MsigDB edges added after merge:     {len(edges_df) - before:,}")
    print(f"Total edges now (unique triples):       {len(edges_df):,}\n")

    # 5) PPI (optional preprocessing)
    print("=== Adding PPI edges (high_confidence_score) via Ensembl -> gene mapping ===")
    ppi_processed = preprocess_ppi(PPI_FILE, min_confidence_score=0)

    ppi_edges = build_ppi_edges(
        ppi_processed,
        GENE_PROTEIN_MAP,
        existing_protein_ids=protein_ids,
        min_confidence_score=0,
    )
    before = len(edges_df)
    edges_df = _merge_edges(edges_df, ppi_edges)
    print(f"PPI produced edges:                     {len(ppi_edges):,}")
    print(f"New PPI edges added after merge:        {len(edges_df) - before:,}")
    print(f"Total edges now (unique triples):       {len(edges_df):,}\n")

    # 6) PTMcode2 (preprocess once, then run builder)
    print("=== Adding PTMcode2 site-site coevolution edges (restricted to existing sites) ===")
    ptmcode_processed = preprocess_ptmcode2(PTMCODE_FILE)

    ptmcode_edges = build_ptmcode2_edges(
        ptmcode_processed,
        GENE_PROTEIN_MAP,
        existing_site_ids=site_ids,
    )
    before = len(edges_df)
    edges_df = _merge_edges(edges_df, ptmcode_edges)
    print(f"PTMcode2 produced edges:                {len(ptmcode_edges):,}")
    print(f"New PTMcode2 edges added after merge:   {len(edges_df) - before:,}")
    print(f"Total edges now (unique triples):       {len(edges_df):,}\n")

    # Write outputs
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    nodes_out = os.path.join(OUTPUT_DIR, "nodes.csv.gz")
    edges_out = os.path.join(OUTPUT_DIR, "edges.csv.gz")

    nodes_df.to_csv(nodes_out, index=False, compression="gzip")
    edges_df.to_csv(edges_out, index=False, compression="gzip")

    print("Wrote outputs:")
    print(f"  {nodes_out}")
    print(f"  {edges_out}")
    print("✅ Done.")


if __name__ == "__main__":
    main()
