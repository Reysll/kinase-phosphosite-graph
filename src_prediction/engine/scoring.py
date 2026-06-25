from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Dict, List

import numpy as np
import pandas as pd


def _cosine_similarity(a: np.ndarray, b: np.ndarray) -> float:
    denom = np.linalg.norm(a) * np.linalg.norm(b)
    if denom == 0:
        return float("nan")
    return float(np.dot(a, b) / denom)


class PairScorer(ABC):
    """
    Abstract scorer interface.

    Later, a supervised model can implement this same interface.
    """

    @abstractmethod
    def score_site_against_kinases(
        self,
        site_node_id: str,
        candidate_kinases: List[str],
        embeddings: Dict[str, np.ndarray],
    ) -> pd.DataFrame:
        """
        Return a DataFrame with columns:
          - kinase_node_id
          - site_node_id
          - score
        """
        raise NotImplementedError


class CosineSimilarityScorer(PairScorer):
    """
    Baseline direct scorer:
    score(kinase, site) = cosine(embedding_kinase, embedding_site)
    """

    def score_site_against_kinases(
        self,
        site_node_id: str,
        candidate_kinases: List[str],
        embeddings: Dict[str, np.ndarray],
    ) -> pd.DataFrame:
        if site_node_id not in embeddings:
            return pd.DataFrame(columns=["kinase_node_id", "site_node_id", "score"])

        site_vec = embeddings[site_node_id]
        rows = []

        for kinase_id in candidate_kinases:
            if kinase_id not in embeddings:
                continue
            score = _cosine_similarity(embeddings[kinase_id], site_vec)
            rows.append(
                {
                    "kinase_node_id": kinase_id,
                    "site_node_id": site_node_id,
                    "score": score,
                }
            )

        out = pd.DataFrame(rows)
        if out.empty:
            return out

        out = out.sort_values(
            by=["score", "kinase_node_id"],
            ascending=[False, True],
        ).reset_index(drop=True)

        return out