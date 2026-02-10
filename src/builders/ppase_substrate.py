# src/builders/ppase_substrate.py
from __future__ import annotations

import re
from typing import List, Tuple

import pandas as pd

from src.io_utils import read_table
from src.normalize import normalize_protein, normalize_site, make_site_id


# Regex to detect site tokens (Ser/Thr/Tyr or S/T/Y + number)
_SITE_PATTERN = re.compile(r"(SER|THR|TYR|S|T|Y)[\s\-]*([0-9]+)", re.IGNORECASE)


def _has_exactly_one_site_token(site_raw: str) -> bool:
    """
    Keep only rows where exactly one site is reported.
    Rows with multiple sites are dropped intentionally.
    """
    if not site_raw or not isinstance(site_raw, str):
        return False
    matches = _SITE_PATTERN.findall(site_raw)
    return len(matches) == 1


def build_ppase_substrate_graph(path: str) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build nodes and edges from PPase_protSubtrates_201903.xls.

    Columns (1-indexed from spec):
      - B: phosphatase
      - D: substrate
      - E: site (Ser/Thr/Tyr naming)

    Nodes:
      - phosphatase (node_type = 'phosphatase')
      - protein (node_type = 'protein')
      - site (node_type = 'site', id = PROTEIN-S88)

    Edges:
      - protein -> site : 'has_site'
      - phosphatase -> site : 'dephosphorylates'
    """

    df = read_table(path)

    # Column indices (0-based)
    phosph_col = 1   # B
    substrate_col = 3  # D
    site_col = 4    # E

    nodes: List[Tuple[str, str]] = []
    edges: List[Tuple[str, str, str]] = []

    # Counters for transparency
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

        phosph_str = str(phosph_raw).strip()
        sub_str = str(sub_raw).strip()
        site_str = str(site_raw).strip()

        # Drop rows with multiple sites
        if not _has_exactly_one_site_token(site_str):
            rows_dropped_multi_site += 1
            continue

        phosph = normalize_protein(phosph_str)
        substrate = normalize_protein(sub_str)
        site_norm = normalize_site(site_str)

        if not phosph or not substrate:
            rows_dropped_protein_parse += 1
            continue

        if site_norm is None:
            rows_dropped_site_parse += 1
            continue

        site_id = make_site_id(substrate, site_norm)

        # Nodes
        nodes.append((phosph, "phosphatase"))
        nodes.append((substrate, "protein"))
        nodes.append((site_id, "site"))

        # Edges
        edges.append((substrate, site_id, "has_site"))
        edges.append((phosph, site_id, "dephosphorylates"))

        rows_kept += 1

    nodes_df = pd.DataFrame(nodes, columns=["node_id", "node_type"]).drop_duplicates()
    edges_df = pd.DataFrame(edges, columns=["source", "target", "relation"]).drop_duplicates()

    # Final stats print
    print("PPase parsing stats:")
    print(f"  rows_seen={rows_seen:,}")
    print(f"  rows_kept={rows_kept:,}")
    print(f"  rows_dropped_multi_site={rows_dropped_multi_site:,}")
    print(f"  rows_dropped_missing_values={rows_dropped_missing_values:,}")
    print(f"  rows_dropped_site_parse={rows_dropped_site_parse:,}")
    print(f"  rows_dropped_protein_parse={rows_dropped_protein_parse:,}")

    return nodes_df, edges_df
