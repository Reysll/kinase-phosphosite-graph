# src/main.py
from __future__ import annotations

import os
import pandas as pd

from src.builders.kinase_substrate import build_kinase_substrate_graph
from src.builders.ppase_substrate import build_ppase_substrate_graph
from src.builders.ptmsigdb import build_ptmsigdb_edges
from src.sanity_checks import check_site_collisions
from src.config import KINASE_SUBSTRATE_FILE, PPASE_FILE, PTMSIGDB_FILE, OUTPUT_DIR


def merge_nodes(nodes_a: pd.DataFrame, nodes_b: pd.DataFrame) -> pd.DataFrame:
    """
    Merge node tables and resolve node_type collisions.

    Rule:
    - If a node_id appears with multiple node_type values,
      coerce node_type to 'enzyme'.
    """
    merged = pd.concat([nodes_a, nodes_b], ignore_index=True)

    type_counts = merged.groupby("node_id")["node_type"].nunique()
    collided_ids = set(type_counts[type_counts > 1].index)

    if collided_ids:
        print(f"WARNING: node_type collisions detected: {len(collided_ids):,} node_ids")
        print("Example node_type collisions:", list(sorted(collided_ids))[:10])
        merged.loc[merged["node_id"].isin(collided_ids), "node_type"] = "enzyme"

    return merged.drop_duplicates(subset=["node_id"], keep="first")




def merge_edges(edges_a: pd.DataFrame, edges_b: pd.DataFrame) -> pd.DataFrame:
    return (
        pd.concat([edges_a, edges_b], ignore_index=True)
        .drop_duplicates(subset=["source", "target", "relation"])
    )


def main() -> None:
    # 1) Base: kinase-substrate graph
    nodes_df, edges_df = build_kinase_substrate_graph(KINASE_SUBSTRATE_FILE)

    # 2) Add: PPase substrate graph
    pp_nodes, pp_edges = build_ppase_substrate_graph(PPASE_FILE)

    nodes_df = merge_nodes(nodes_df, pp_nodes)
    edges_df = merge_edges(edges_df, pp_edges)

    # 3) Sanity checks on merged graph
    check_site_collisions(nodes_df)

    # 4) PTMsigDB site-site edges only for existing sites
    site_ids = set(nodes_df.loc[nodes_df["node_type"] == "site", "node_id"])
    ptm_edges = build_ptmsigdb_edges(PTMSIGDB_FILE, site_ids)
    print(f"PTMsigDB edges added: {len(ptm_edges):,}")
    edges_df = merge_edges(edges_df, ptm_edges)

    # 5) Write outputs
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    nodes_df.to_csv(os.path.join(OUTPUT_DIR, "nodes.csv.gz"), index=False, compression="gzip")
    edges_df.to_csv(os.path.join(OUTPUT_DIR, "edges.csv.gz"), index=False, compression="gzip")

    print("Wrote outputs:")
    print(" ", os.path.join(OUTPUT_DIR, "nodes.csv.gz"))
    print(" ", os.path.join(OUTPUT_DIR, "edges.csv.gz"))
    print("✅ Done.")


if __name__ == "__main__":
    main()
