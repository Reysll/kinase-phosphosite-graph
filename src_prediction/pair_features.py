from __future__ import annotations

from typing import Dict, Iterable, List

import numpy as np
import pandas as pd


def build_pair_feature_vector(kinase_vec: np.ndarray, site_vec: np.ndarray) -> np.ndarray:
    """
    Pair features:
      [kinase_vec, site_vec, abs(kinase-site), kinase*site]
    """
    return np.concatenate(
        [
            kinase_vec,
            site_vec,
            np.abs(kinase_vec - site_vec),
            kinase_vec * site_vec,
        ]
    )


def build_pair_feature_table(
    pairs_df: pd.DataFrame,
    embeddings: Dict[str, np.ndarray],
    kinase_col: str = "kinase_node_id",
    site_col: str = "site_node_id",
    label_col: str | None = None,
) -> pd.DataFrame:
    rows = []

    for row in pairs_df.itertuples(index=False):
        kinase_id = str(getattr(row, kinase_col))
        site_id = str(getattr(row, site_col))

        if kinase_id not in embeddings or site_id not in embeddings:
            continue

        feat = build_pair_feature_vector(embeddings[kinase_id], embeddings[site_id])

        out_row = {
            "kinase_node_id": kinase_id,
            "site_node_id": site_id,
        }

        for i, val in enumerate(feat):
            out_row[f"f_{i}"] = float(val)

        if label_col is not None:
            out_row[label_col] = int(getattr(row, label_col))

        rows.append(out_row)

    return pd.DataFrame(rows)