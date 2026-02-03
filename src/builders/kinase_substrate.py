# src/builders/kinase_substrate.py
from typing import Dict, List, Set, Tuple

import pandas as pd

from src.normalize import normalize_protein, normalize_site, make_site_id
from src.io_utils import read_table



def build_kinase_substrate_graph(
    path: str,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build nodes and edges from the Kinase_Substrate_Dataset.

    Returns:
        nodes_df: node_id, node_type
        edges_df: source, target, relation
    """

    df = read_table(path)

    # Column indices based on dataset description
    kinase_col = df.columns[0]     # A
    kinase_org_col = df.columns[3] # D
    substrate_col = df.columns[7]  # H
    substrate_org_col = df.columns[8]  # I
    site_col = df.columns[9]       # J

    # Filter human–human
    human = (
        df[kinase_org_col].str.lower().str.contains("human")
        & df[substrate_org_col].str.lower().str.contains("human")
    )
    df = df[human]

    nodes: Dict[str, str] = {}
    edges: List[Tuple[str, str, str]] = []

    for _, row in df.iterrows():
        kinase = normalize_protein(row[kinase_col])
        protein = normalize_protein(row[substrate_col])

        site = normalize_site(row[site_col])
        if site is None:
            continue

        site_id = make_site_id(protein, site)

        nodes[kinase] = "kinase"
        nodes[protein] = "protein"
        nodes[site_id] = "site"

        edges.append((protein, site_id, "has_site"))
        edges.append((kinase, site_id, "phosphorylates"))

    nodes_df = pd.DataFrame(
        [{"node_id": k, "node_type": v} for k, v in nodes.items()]
    )

    edges_df = pd.DataFrame(
        edges, columns=["source", "target", "relation"]
    )

    return nodes_df, edges_df
