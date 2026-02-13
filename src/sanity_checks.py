# src/sanity_checks.py
from __future__ import annotations

import pandas as pd

from src.ids import parse_site_node_id


def check_site_collisions(nodes_df: pd.DataFrame) -> None:
    """
    Sanity check for SITE nodes.

    Under Option B, site node IDs look like:
      SITE:MAP2K1-S298

    We want to confirm:
      1) Many site labels (e.g., S12) appear across multiple proteins (expected)
      2) No duplicate (protein, site_label) within the same protein
    """
    sites = nodes_df[nodes_df["node_type"] == "site"].copy()

    proteins: list[str] = []
    labels: list[str] = []
    bad_parse = 0

    for node_id in sites["node_id"].tolist():
        parsed = parse_site_node_id(str(node_id))
        if parsed is None:
            bad_parse += 1
            proteins.append("PARSE_FAIL")
            labels.append("PARSE_FAIL")
            continue
        proteins.append(parsed.gene)
        labels.append(parsed.site_label)

    sites["protein"] = proteins
    sites["site_label"] = labels

    if bad_parse > 0:
        print(f"WARNING: {bad_parse:,} site node_ids could not be parsed as SITE:<PROTEIN>-<LABEL>")

    # How many proteins share a given site label (expected to be >0 for many labels)
    shared_counts = sites.groupby("site_label")["protein"].nunique()
    shared_labels = shared_counts[shared_counts > 1].sort_values(ascending=False)

    print(f"Total site nodes: {len(sites):,}")
    print(f"Site labels shared across proteins (expected): {len(shared_labels):,}")

    if len(shared_labels) > 0:
        print("Examples (site_label -> #proteins):")
        for label, nprot in shared_labels.head(10).items():
            print(f"  {label} -> {nprot}")

    # Duplicate within-protein sites should not exist
    dup_within = sites.groupby(["protein", "site_label"]).size()
    dup_within = dup_within[dup_within > 1]

    if not dup_within.empty:
        print("WARNING: duplicate site nodes within same protein detected!")
        print(dup_within.head(10))
    else:
        print("No duplicate sites within proteins.")
