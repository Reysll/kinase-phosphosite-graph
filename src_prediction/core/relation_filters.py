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

# PhosphoSitePlus KSA edges — these are the prediction TARGETS and must be
# excluded from the graph used to compute node2vec embeddings (leakage fix).
PSP_KSA_RELATIONS: Set[str] = {"phosphorylates"}

# Relations kept in the embedding graph (everything except PSP KSA edges)
EMBEDDING_RELATIONS: Set[str] = GENERIC_BASE_RELATIONS - PSP_KSA_RELATIONS

# --- Liver context correlation edge sets ---
# FC-based (fold-change across tumor samples vs control mean) — original
SITE_CORR_FC_RELATIONS: Set[str] = {"site_corr_fc_pos", "site_corr_fc_neg"}

# Control-specific (raw abundance correlations across healthy/control samples)
SITE_CORR_CTRL_RELATIONS: Set[str] = {"site_corr_ctrl_pos", "site_corr_ctrl_neg"}

# Cancer-specific (raw abundance correlations across tumor/cancer samples)
SITE_CORR_CANCER_RELATIONS: Set[str] = {"site_corr_cancer_pos", "site_corr_cancer_neg"}

# Legacy alias — kept for backward compatibility with existing scripts
SITE_CORR_RELATIONS: Set[str] = (
    SITE_CORR_FC_RELATIONS | SITE_CORR_CTRL_RELATIONS | SITE_CORR_CANCER_RELATIONS
)

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
