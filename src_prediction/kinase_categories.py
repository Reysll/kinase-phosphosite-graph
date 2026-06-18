"""
Substrate-count-based kinase category assignment for within-category ranking.

Motivation (Dr. Ayati, 2026-06-17):
  When all 420 candidate kinases compete in a single ranked list, "poor"
  kinases (few known substrates) cannot realistically outrank "rich" ones
  (hundreds of substrates), because the model is trained on more positive
  examples for rich kinases. Category-based ranking removes that asymmetry:
  each site gets one top-1 prediction per category.

Design decisions vs. the supervisor's request
----------------------------------------------
1. Pre-filter counts (our choice):
   Dr. Ayati left open "before filtering for multi-kinase or after". We use
   counts from the FULL generic graph (all PSP edges, not just multi-kinase
   sites). A kinase's richness should reflect its global PSP annotation;
   using post-filter counts would make a kinase "poor" simply because most of
   its substrates are single-kinase sites (excluded from our evaluation).

2. Fixed thresholds (per Dr. Ayati's suggestion):
   POOR_MAX=5 and RICH_MIN=20 were suggested by Dr. Ayati as a starting point.
   With our 420 candidates they give a near-even split: 166 poor / 133 average /
   121 rich. Percentile-based thresholds are a natural next step if the
   supervisor wants strictly equal-sized bins.

3. Category hub caveat:
   The known hub kinases land in different categories — VRK3 (1 substrate) is
   poor, VRK2 (15) is average, ULK1 (77) is rich. Category-based ranking
   therefore reduces CROSS-category competition but does NOT eliminate
   within-category hub bias. SpectralEmbeddingStrategy (Option 1) is needed
   to address the embedding-level bias within each category.
"""

from __future__ import annotations

from typing import Dict, List

import pandas as pd

# Provisional thresholds (Dr. Ayati, 2026-06-17). Produces 166/133/121 split
# on our 420-kinase candidate set. Revisit if the supervisor prefers equal bins.
POOR_MAX: int = 5
RICH_MIN: int = 20

CATEGORIES = ("poor", "average", "rich")


def compute_psp_substrate_counts(generic_edges: pd.DataFrame) -> pd.Series:
    """
    Count unique PSP substrate sites per kinase from the full generic edge table.

    Returns pd.Series indexed by kinase node_id (PROTEIN:GENENAME), sorted
    descending. Kinases with no 'phosphorylates' edges are absent from the result
    and will be assigned to 'poor' by assign_kinase_categories.
    """
    psp = generic_edges.loc[generic_edges["relation"] == "phosphorylates"]
    return psp.groupby("source")["target"].nunique().sort_values(ascending=False)


def assign_kinase_categories(
    kinase_node_ids: List[str],
    substrate_counts: pd.Series,
) -> Dict[str, str]:
    """
    Map each kinase to 'poor', 'average', or 'rich'.

    Kinases absent from substrate_counts default to 'poor'.
    """
    counts = substrate_counts.to_dict()
    result: Dict[str, str] = {}
    for kid in kinase_node_ids:
        n = counts.get(kid, 0)
        if n <= POOR_MAX:
            result[kid] = "poor"
        elif n <= RICH_MIN:
            result[kid] = "average"
        else:
            result[kid] = "rich"
    return result


def print_category_summary(
    kinase_categories: Dict[str, str],
    substrate_counts: pd.Series,
) -> None:
    """Print category sizes and known hub-kinase placement for quick sanity check."""
    sizes = {cat: sum(1 for c in kinase_categories.values() if c == cat) for cat in CATEGORIES}
    total = sum(sizes.values())
    print("Kinase category distribution (thresholds: poor <=5, average 5-20, rich >20):")
    for cat in CATEGORIES:
        print(f"  {cat:8s}: {sizes[cat]:3d} kinases ({sizes[cat] / total * 100:.1f}%)")
    print()
    # Show placement of known hub kinases
    hub_kinases = ["PROTEIN:ULK1", "PROTEIN:VRK2", "PROTEIN:VRK3", "PROTEIN:WNK1"]
    print("Hub kinase categories (context for interpreting top-1 predictions):")
    counts_dict = substrate_counts.to_dict()
    for k in hub_kinases:
        cat = kinase_categories.get(k, "not a candidate")
        n = counts_dict.get(k, 0)
        print(f"  {k}: {n} substrates -> {cat}")
    print()
