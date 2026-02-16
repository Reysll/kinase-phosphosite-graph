# src/builders/ppase_substrate.py
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


def build_ppase_substrate_graph(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Unified protein-level design.

    Nodes:
      - PROTEIN:<phosphatase>
      - PROTEIN:<substrate>
      - SITE:<substrate>-<site>

    Edges:
      - PROTEIN:<substrate>    -> SITE   relation="has_site"
      - PROTEIN:<phosphatase>  -> SITE   relation="dephosphorylates"
    """
    df = read_table(path)

    phosph_col = 1
    substrate_col = 3
    site_col = 4

    nodes: List[Tuple[str, str]] = []
    edges: List[Tuple[str, str, str]] = []

    rows_seen = 0
    rows_kept = 0
    rows_dropped_multi_site = 0
    rows_dropped_missing_values = 0
    rows_dropped_site_parse = 0
    rows_dropped_protein_parse = 0

    for _, row in df.iterrows():
        rows_seen += 1

        try:
            phosph_raw = row.iloc[phosph_col]
            sub_raw = row.iloc[substrate_col]
            site_raw = row.iloc[site_col]
        except IndexError:
            rows_dropped_missing_values += 1
            continue

        if pd.isna(phosph_raw) or pd.isna(sub_raw) or pd.isna(site_raw):
            rows_dropped_missing_values += 1
            continue

        site_str = str(site_raw).strip()
        if not _has_exactly_one_site_token(site_str):
            rows_dropped_multi_site += 1
            continue

        ppase_name = normalize_protein(str(phosph_raw))
        substrate_name = normalize_protein(str(sub_raw))
        site_label = normalize_site(site_str)

        if not ppase_name or not substrate_name:
            rows_dropped_protein_parse += 1
            continue
        if site_label is None:
            rows_dropped_site_parse += 1
            continue

        pp_node = protein_id(ppase_name)
        p_node = protein_id(substrate_name)
        s_node = site_id(substrate_name, site_label)

        nodes.append((pp_node, "protein"))
        nodes.append((p_node, "protein"))
        nodes.append((s_node, "site"))

        edges.append((p_node, s_node, "has_site"))
        edges.append((pp_node, s_node, "dephosphorylates"))

        rows_kept += 1

    nodes_df = pd.DataFrame(nodes, columns=["node_id", "node_type"]).drop_duplicates()
    edges_df = pd.DataFrame(edges, columns=["source", "target", "relation"]).drop_duplicates()

    print("PPase parsing stats:")
    print(f"  rows_seen={rows_seen:,}")
    print(f"  rows_kept={rows_kept:,}")
    print(f"  rows_dropped_multi_site={rows_dropped_multi_site:,}")
    print(f"  rows_dropped_missing_values={rows_dropped_missing_values:,}")
    print(f"  rows_dropped_site_parse={rows_dropped_site_parse:,}")
    print(f"  rows_dropped_protein_parse={rows_dropped_protein_parse:,}")

    return nodes_df, edges_df
