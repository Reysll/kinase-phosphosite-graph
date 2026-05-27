"""
LOO evaluation on the cancer-specific (tumor sample) network.

Correlation edges: site_corr_cancer_pos / site_corr_cancer_neg
These are built from Pearson correlations of raw phosphosite abundances
across the 18 liver cancer / tumor samples in the liver proteomics dataset.

Compare with run_multi_kinase_liver.py (FC-based): both use tumor samples,
but the FC model correlates fold-change vectors (tumor / mean_control) while
this model correlates raw tumor abundances directly.

Run after rebuilding the liver graph (src.run_liver now includes all 3
correlation types in a single Liver_network_edges.csv.gz file).
"""

from __future__ import annotations

from src_prediction.config import (
    CANDIDATE_KINASES_OUT,
    LIVER_EDGES,
    LIVER_NODES,
    PRED_OUTPUTS_DIR,
)
from src_prediction.embedding_strategy import Node2VecStrategy
from src_prediction.experiment_utils import write_experiment_outputs
from src_prediction.graph_loader import load_graph
from src_prediction.io_utils import read_csv_auto
from src_prediction.leave_one_out import run_leave_one_out
from src_prediction.metrics import summarize_results, summarize_results_text
from src_prediction.relation_filters import GENERIC_BASE_RELATIONS, SITE_CORR_CANCER_RELATIONS


def main() -> None:
    experiment_name = "cancer_multi_kinase"

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

    print("=== Loading liver graph (cancer-specific correlation variant) ===")
    graph = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
    print(f"Graph nodes: {len(graph.nodes):,}")
    print(f"Graph edges: {len(graph.edges):,}\n")

    print("=== Loading evaluation inputs ===")
    candidate_kinases = read_csv_auto(CANDIDATE_KINASES_OUT)
    positive_edges = read_csv_auto(frozen_trials_path)
    print(f"Candidate kinases: {len(candidate_kinases):,}")
    print(f"Multi-kinase LOO trials: {len(positive_edges):,}\n")

    # Cancer model: base relations + tumor-sample correlation edges only.
    # PSP 'phosphorylates' edges are automatically excluded from the embedding
    # graph inside run_leave_one_out (leakage fix); they remain as training labels.
    allowed_relations = set(GENERIC_BASE_RELATIONS) | SITE_CORR_CANCER_RELATIONS

    print("=== Running leave-one-out on multi-kinase sites: cancer (tumor) model ===")

    results = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=positive_edges,
        candidate_kinases=candidate_kinases,
        embedding_strategy=strategy,
        allowed_relations=allowed_relations,
        max_trials=None,
        random_state=42,
        verbose_every=100,
        n_jobs_outer=8,
        max_negatives_per_site=50,
    )

    metrics_df = summarize_results(results)
    summary_text = summarize_results_text(results)

    paths = write_experiment_outputs(
        experiment_name=experiment_name,
        base_dir=PRED_OUTPUTS_DIR,
        results_df=results,
        metrics_df=metrics_df,
        summary_text=summary_text,
    )

    print(summary_text)
    print()
    print(f"Wrote: {paths['results']}")
    print(f"Wrote: {paths['metrics']}")
    print(f"Wrote: {paths['summary']}")
    print("Done.")


if __name__ == "__main__":
    main()
