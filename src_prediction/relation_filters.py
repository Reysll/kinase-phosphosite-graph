from __future__ import annotations

from typing import Optional, Set

import pandas as pd


GENERIC_BASE_RELATIONS: Set[str] = {
    "has_site",
    "phosphorylates",
    "dephosphorylates",
    "site_same_pathway",
    "protein_same_pathway",
    "ppi_high_confidence",
    "site_coevolution",
}

SITE_CORR_RELATIONS: Set[str] = {
    "site_corr_fc_pos",
    "site_corr_fc_neg",
}

PROTEIN_CORR_RELATIONS: Set[str] = {
    "protein_corr_fc_pos",
    "protein_corr_fc_neg",
}


def filter_edges_by_relation(
    edges_df: pd.DataFrame,
    allowed_relations: Optional[Set[str]],
) -> pd.DataFrame:
    if allowed_relations is None:
        return edges_df.copy()
    return edges_df.loc[edges_df["relation"].isin(allowed_relations)].copy()