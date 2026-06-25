"""
50-trial smoke test: SpectralEmbeddingStrategy + category-based ranking
(Option 1 + Option 2 combined).

Purpose: verify that the degree-normalized embedding runs correctly and
check whether hub bias is reduced before committing to the full 6,909-trial
run. Expected runtime: ~2-4 min for the spectral decomposition on the liver
graph (~16 K nodes), then ~1 min for 50 trials.

What to look for in the output:
  - Top-1 per category: is it still one kinase for every trial? If spectral
    embedding works, the dominance of ULK1/VRK2/VRK3 should be reduced.
  - adjusted_held_out_rank_in_category vs. node2vec equivalent: lower is
    better; compare with run_test_categorized.py output.
  - If top-1 is still a single kinase for 100% of 50 trials, the hub bias
    survived the embedding swap and we need to investigate graph structure
    (e.g., very high-degree nodes may still dominate via Laplacian structure).

Outputs: outputs_prediction/test_spectral_categorized/
"""

from __future__ import annotations

from collections import Counter

import pandas as pd

from src_prediction.core.config import (
    CANDIDATE_KINASES_OUT,
    GENERIC_EDGES,
    LIVER_EDGES,
    LIVER_NODES,
    PRED_OUTPUTS_DIR,
)
from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy
from src_prediction.core.experiment_utils import write_experiment_outputs
from src_prediction.core.graph_loader import load_graph
from src_prediction.core.io_utils import read_csv_auto
from src_prediction.data_prep.kinase_categories import (
    assign_kinase_categories,
    compute_psp_substrate_counts,
    print_category_summary,
)
from src_prediction.engine.leave_one_out import run_leave_one_out
from src_prediction.core.metrics import summarize_results, summarize_results_text
from src_prediction.core.relation_filters import GENERIC_BASE_RELATIONS, SITE_CORR_FC_RELATIONS

MAX_TRIALS = 50  # increase to None once validated


def _print_category_rank_summary(results: pd.DataFrame) -> None:
    cat_col = "held_out_kinase_category"
    rank_col = "adjusted_held_out_rank_in_category"

    print("=== Within-category adjusted rank ===")
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
    experiment_name = "test_spectral_categorized"

    frozen_trials_path = PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz"

    # SpectralEmbeddingStrategy: degree-normalized via normalized Laplacian.
    # directed=False is required — spectral embedding needs a symmetric adjacency.
    # The LOO runner reads strategy.directed before building the embedding graph,
    # so no manual override is needed in allowed_relations.
    strategy = SpectralEmbeddingStrategy(
        dimensions=32,
        directed=False,
        random_state=42,
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

    kinase_node_ids = candidate_kinases["node_id"].astype(str).tolist()
    kinase_categories = assign_kinase_categories(kinase_node_ids, substrate_counts)
    print_category_summary(kinase_categories, substrate_counts)

    # Spectral embedding note: PSP 'phosphorylates' edges are still excluded
    # from the embedding graph inside run_leave_one_out (leakage fix applies
    # regardless of embedding strategy).
    allowed_relations = set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS

    print("=== Running LOO (spectral + category ranking, liver FC, 50 trials) ===")
    print("    Note: spectral decomposition on ~16 K nodes takes ~2-4 min.\n")
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
