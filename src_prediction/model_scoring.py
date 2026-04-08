from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from src_prediction.pair_features import build_pair_feature_table
from src_prediction.negative_sampling import sample_negative_pairs_for_site


@dataclass
class TrainedPairModelResult:
    model: LogisticRegression
    train_table: pd.DataFrame


def train_logistic_pair_model(
    train_positive_edges: pd.DataFrame,
    candidate_kinase_ids: List[str],
    site_to_true_kinases: Dict[str, List[str]],
    embeddings: Dict[str, np.ndarray],
    max_negatives_per_site: int = 50,
) -> TrainedPairModelResult:
    """
    Train a simple logistic regression model for kinase-site scoring.

    Positives:
      known kinase-site edges in training folds

    Negatives:
      for each positive site, candidate kinases not known to be true for that site
    """
    pos = train_positive_edges.copy()
    pos["label"] = 1

    neg_tables = []
    seen_sites = sorted(set(pos["site_node_id"].astype(str)))

    for site in seen_sites:
        true_set = set(site_to_true_kinases.get(site, []))
        neg_df = sample_negative_pairs_for_site(
            site_node_id=site,
            candidate_kinase_ids=candidate_kinase_ids,
            true_kinases_for_site=true_set,
            max_negatives=max_negatives_per_site,
        )
        neg_tables.append(neg_df)

    neg = pd.concat(neg_tables, ignore_index=True) if neg_tables else pd.DataFrame(
        columns=["kinase_node_id", "site_node_id", "label"]
    )

    train_pairs = pd.concat(
        [pos[["kinase_node_id", "site_node_id", "label"]], neg],
        ignore_index=True,
    )

    feature_table = build_pair_feature_table(
        pairs_df=train_pairs,
        embeddings=embeddings,
        label_col="label",
    )

    feature_cols = [c for c in feature_table.columns if c.startswith("f_")]
    X = feature_table[feature_cols].values
    y = feature_table["label"].values

    model = LogisticRegression(
        max_iter=1000,
        class_weight="balanced",
        solver="liblinear",
        random_state=42,
    )
    model.fit(X, y)

    return TrainedPairModelResult(model=model, train_table=feature_table)


def score_site_with_trained_model(
    model: LogisticRegression,
    site_node_id: str,
    candidate_kinase_ids: List[str],
    embeddings: Dict[str, np.ndarray],
) -> pd.DataFrame:
    pairs = pd.DataFrame(
        {
            "kinase_node_id": candidate_kinase_ids,
            "site_node_id": [site_node_id] * len(candidate_kinase_ids),
        }
    )

    feature_table = build_pair_feature_table(
        pairs_df=pairs,
        embeddings=embeddings,
        label_col=None,
    )

    if feature_table.empty:
        return pd.DataFrame(columns=["kinase_node_id", "site_node_id", "score"])

    feature_cols = [c for c in feature_table.columns if c.startswith("f_")]
    X = feature_table[feature_cols].values
    probs = model.predict_proba(X)[:, 1]

    out = feature_table[["kinase_node_id", "site_node_id"]].copy()
    out["score"] = probs
    out = out.sort_values(
        by=["score", "kinase_node_id"],
        ascending=[False, True],
    ).reset_index(drop=True)

    return out