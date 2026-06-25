from __future__ import annotations

from typing import Dict, List, Set

import pandas as pd


def sample_negative_pairs_for_site(
    site_node_id: str,
    candidate_kinase_ids: List[str],
    true_kinases_for_site: Set[str],
    max_negatives: int | None = None,
) -> pd.DataFrame:
    """
    Build negative kinase-site pairs for one site by excluding known true kinases.
    """
    negatives = [k for k in candidate_kinase_ids if k not in true_kinases_for_site]

    if max_negatives is not None:
        negatives = negatives[:max_negatives]

    return pd.DataFrame(
        {
            "kinase_node_id": negatives,
            "site_node_id": [site_node_id] * len(negatives),
            "label": [0] * len(negatives),
        }
    )