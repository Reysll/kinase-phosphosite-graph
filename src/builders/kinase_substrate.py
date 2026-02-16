# src/builders/kinase_substrate.py
from __future__ import annotations

import re
from typing import List, Tuple

import pandas as pd

from src.io_utils import read_table
from src.normalize import normalize_protein, normalize_site
from src.ids import protein_id, site_id


_SITE_PATTERN = re.compile(r"(SER|THR|TYR|S|T|Y)[\s\-]*([0-9]+)", re.IGNORECASE)


def _has_exactly_one_site_token(site_raw: str) -> bool:
    if not site_raw or not isinstance(site_raw, str):
        return False
    matches = _SITE_PATTERN.findall(site_raw)
    return len(matches) == 1


def build_kinase_substrate_graph(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Unified protein-level design.

    Nodes:
      - PROTEIN:<kinase>
      - PROTEIN:<substrate>
      - SITE:<substrate>-<site>

    Edges:
      - PROTEIN:<substrate> -> SITE   relation="has_site"
      - PROTEIN:<kinase>    -> SITE   relation="phosphorylates"
    """
    df = read_table(path)

    kinase_col = 0
    kinase_org_col = 3
    substrate_col = 7
    substrate_org_col = 8
    site_col = 9

    nodes: List[Tuple[str, str]] = []
    edges: List[Tuple[str, str, str]] = []

    rows_seen = 0
    rows_kept = 0
    rows_dropped_nonhuman = 0
    rows_dropped_missing_values = 0
    rows_dropped_multi_site = 0
    rows_dropped_parse = 0

    for _, row in df.iterrows():
        rows_seen += 1

        try:
            kinase_raw = row.iloc[kinase_col]
            kinase_org_raw = row.iloc[kinase_org_col]
            substrate_raw = row.iloc[substrate_col]
            substrate_org_raw = row.iloc[substrate_org_col]
            site_raw = row.iloc[site_col]
        except IndexError:
            rows_dropped_missing_values += 1
            continue

        if (
            pd.isna(kinase_raw)
            or pd.isna(kinase_org_raw)
            or pd.isna(substrate_raw)
            or pd.isna(substrate_org_raw)
            or pd.isna(site_raw)
        ):
            rows_dropped_missing_values += 1
            continue

        kinase_org = str(kinase_org_raw).strip().lower()
        substrate_org = str(substrate_org_raw).strip().lower()

        if "human" not in kinase_org or "human" not in substrate_org:
            rows_dropped_nonhuman += 1
            continue

        site_str = str(site_raw).strip()
        if not _has_exactly_one_site_token(site_str):
            rows_dropped_multi_site += 1
            continue

        kinase_name = normalize_protein(str(kinase_raw))
        substrate_name = normalize_protein(str(substrate_raw))
        site_label = normalize_site(site_str)

        if not kinase_name or not substrate_name or site_label is None:
            rows_dropped_parse += 1
            continue

        k_node = protein_id(kinase_name)
        p_node = protein_id(substrate_name)
        s_node = site_id(substrate_name, site_label)

        nodes.append((k_node, "protein"))
        nodes.append((p_node, "protein"))
        nodes.append((s_node, "site"))

        edges.append((p_node, s_node, "has_site"))
        edges.append((k_node, s_node, "phosphorylates"))

        rows_kept += 1

    nodes_df = pd.DataFrame(nodes, columns=["node_id", "node_type"]).drop_duplicates()
    edges_df = pd.DataFrame(edges, columns=["source", "target", "relation"]).drop_duplicates()

    print("Kinase_Substrate parsing stats:")
    print(f"  rows_seen={rows_seen:,}")
    print(f"  rows_kept={rows_kept:,}")
    print(f"  rows_dropped_nonhuman={rows_dropped_nonhuman:,}")
    print(f"  rows_dropped_missing_values={rows_dropped_missing_values:,}")
    print(f"  rows_dropped_multi_site={rows_dropped_multi_site:,}")
    print(f"  rows_dropped_parse={rows_dropped_parse:,}")

    return nodes_df, edges_df