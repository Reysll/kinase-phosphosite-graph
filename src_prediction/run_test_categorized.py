"""
50-trial smoke test: node2vec + category-based ranking (Option 2).

Purpose: verify that within-category ranking is wired correctly before
committing to the full 6,909-trial run. Checks:
  1. Category sizes look right (expect ~166 poor / 133 average / 121 rich).
  2. Top-1 per category shows whether within-category hub bias is still present.
     Known hubs per category: VRK3 (poor, 1 substrate), VRK2 (average, 15),
     ULK1 (rich, 77). Category ranking reduces CROSS-category competition but
     does NOT remove within-category hub bias — that requires Option 1
     (SpectralEmbeddingStrategy, see run_test_spectral.py).
  3. adjusted_held_out_rank_in_category is smaller than adjusted_held_out_rank
     (denominator is ~166/133/121 instead of 420), so raw numbers are not
     directly comparable to previous results.

Outputs: outputs_prediction/test_node2vec_categorized/
"""

from __future__ import annotations

from collections import Counter

import pandas as pd

from src_prediction.config import (
    CANDIDATE_KINASES_OUT,
    GENERIC_EDGES,
    LIVER_EDGES,
    LIVER_NODES,
    PRED_OUTPUTS_DIR,
)
from src_prediction.embedding_strategy import Node2VecStrategy
from src_prediction.experiment_utils import write_experiment_outputs
from src_prediction.graph_loader import load_graph
from src_prediction.io_utils import read_csv_auto
from src_prediction.kinase_categories import (
    assign_kinase_categories,
    compute_psp_substrate_counts,
    print_category_summary,
)
from src_prediction.leave_one_out import run_leave_one_out
from src_prediction.metrics import summarize_results, summarize_results_text
from src_prediction.relation_filters import GENERIC_BASE_RELATIONS, SITE_CORR_FC_RELATIONS

MAX_TRIALS = 50  # increase to 6909 (or None) once validated


def _print_category_rank_summary(results: pd.DataFrame) -> None:
    """Print per-category adjusted rank stats and top-1 frequency."""
    cat_col = "held_out_kinase_category"
    rank_col = "adjusted_held_out_rank_in_category"

    print("=== Within-category adjusted rank (denominator = category size, not 420) ===")
    for cat in ("poor", "average", "rich"):
        subset = results[results[cat_col] == cat][rank_col].dropna()
        if subset.empty:
            print(f"  {cat}: no trials")
            continue
        print(
            f"  {cat:8s}: n={len(subset):3d}  mean={subset.mean():.1f}  median={subset.median():.1f}"
        )
    print()

    print("=== Top-1 predictions per category (% of trials) ===")
    for cat in ("poor", "average", "rich"):
        col = f"top1_{cat}"
        if col not in results.columns:
            continue
        counts = Counter(results[col].dropna())
        total = len(results[col].dropna())
        print(f"  {cat}:")
        for kinase, n in counts.most_common(5):
            print(f"    {kinase}: {n}/{total} ({n / total * 100:.1f}%)")
    print()


def main() -> None:
    experiment_name = "test_node2vec_categorized"

    frozen_trials_path = PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz"

    strategy = Node2VecStrategy(
        dimensions=32,
        walk_length=10,
        num_walks=25,
        workers=4,
        p=1.0,
        q=1.0,
        window=5,
        min_count=1,
        batch_words=10000,
        seed=42,
        directed=True,
    )

    print("=== Loading generic edges for substrate counts ===")
    generic_edges = read_csv_auto(GENERIC_EDGES)
    substrate_counts = compute_psp_substrate_counts(generic_edges)

    print("=== Loading liver graph ===")
    graph = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
    print(f"Graph nodes: {len(graph.nodes):,}")
    print(f"Graph edges: {len(graph.edges):,}\n")

    print("=== Loading evaluation inputs ===")
    candidate_kinases = read_csv_auto(CANDIDATE_KINASES_OUT)
    positive_edges = read_csv_auto(frozen_trials_path)
    print(f"Candidate kinases: {len(candidate_kinases):,}")
    print(f"Multi-kinase LOO trials (full set): {len(positive_edges):,}")
    print(f"This test run: {MAX_TRIALS} trials\n")

    # Build category map from full generic PSP edges (pre-filter counts)
    kinase_node_ids = candidate_kinases["node_id"].astype(str).tolist()
    kinase_categories = assign_kinase_categories(kinase_node_ids, substrate_counts)
    print_category_summary(kinase_categories, substrate_counts)

    allowed_relations = set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS

    print("=== Running LOO (node2vec + category ranking, liver FC, 50 trials) ===")
    results = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=positive_edges,
        candidate_kinases=candidate_kinases,
        embedding_strategy=strategy,
        allowed_relations=allowed_relations,
        max_trials=MAX_TRIALS,
        random_state=42,
        verbose_every=10,
        n_jobs_outer=8,
        max_negatives_per_site=50,
        kinase_categories=kinase_categories,
    )

    print("\n=== Standard metrics (global ranking, 420 candidates) ===")
    summary_text = summarize_results_text(results)
    print(summary_text)
    print()

    _print_category_rank_summary(results)

    metrics_df = summarize_results(results)
    paths = write_experiment_outputs(
        experiment_name=experiment_name,
        base_dir=PRED_OUTPUTS_DIR,
        results_df=results,
        metrics_df=metrics_df,
        summary_text=summary_text,
    )

    print(f"Wrote: {paths['results']}")
    print(f"Wrote: {paths['metrics']}")
    print(f"Wrote: {paths['summary']}")
    print("Done.")


if __name__ == "__main__":
    main()
