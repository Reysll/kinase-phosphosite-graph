from __future__ import annotations

from typing import Dict

import numpy as np
import pandas as pd


def build_pair_feature_table(
    pairs_df: pd.DataFrame,
    embeddings: Dict[str, np.ndarray],
    kinase_col: str = "kinase_node_id",
    site_col: str = "site_node_id",
    label_col: str | None = None,
) -> pd.DataFrame:
    kinase_ids = pairs_df[kinase_col].astype(str).tolist()
    site_ids = pairs_df[site_col].astype(str).tolist()

    valid_indices = [
        i for i, (k, s) in enumerate(zip(kinase_ids, site_ids))
        if k in embeddings and s in embeddings
    ]

    if not valid_indices:
        return pd.DataFrame(columns=["kinase_node_id", "site_node_id"])

    valid_kinases = [kinase_ids[i] for i in valid_indices]
    valid_sites = [site_ids[i] for i in valid_indices]

    K = np.array([embeddings[k] for k in valid_kinases])
    S = np.array([embeddings[s] for s in valid_sites])
    features = np.hstack([K, S, np.abs(K - S), K * S])

    feat_df = pd.DataFrame(features, columns=[f"f_{i}" for i in range(features.shape[1])])
    out = pd.concat(
        [
            pd.DataFrame({"kinase_node_id": valid_kinases, "site_node_id": valid_sites}),
            feat_df,
        ],
        axis=1,
    )

    if label_col is not None:
        out[label_col] = pairs_df[label_col].values[valid_indices]

    return out