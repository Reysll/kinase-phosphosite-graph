from __future__ import annotations

from typing import Dict, List

import pandas as pd


def build_site_to_true_kinases(positive_edges: pd.DataFrame) -> Dict[str, List[str]]:
    """
    Build mapping:
      site_node_id -> sorted list of all known true kinase_node_id values
    using the positive evaluation edge table.

    positive_edges must contain:
      - kinase_node_id
      - site_node_id
    """
    required = {"kinase_node_id", "site_node_id"}
    missing = required - set(positive_edges.columns)
    if missing:
        raise ValueError(f"positive_edges missing required columns: {sorted(missing)}")

    grouped = (
        positive_edges.groupby("site_node_id")["kinase_node_id"]
        .apply(lambda s: sorted(set(map(str, s))))
        .to_dict()
    )
    return grouped