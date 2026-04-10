from __future__ import annotations

from src_prediction.config import (
    GENERIC_EDGES,
    GENERIC_NODES,
    CANDIDATE_KINASES_OUT,
    PRED_OUTPUTS_DIR,
)
from src_prediction.experiment_utils import write_experiment_outputs
from src_prediction.graph_loader import load_graph
from src_prediction.io_utils import read_csv_auto
from src_prediction.leave_one_out import Node2VecParams, run_leave_one_out
from src_prediction.metrics import summarize_results, summarize_results_text
from src_prediction.relation_filters import GENERIC_BASE_RELATIONS


def main() -> None:
    experiment_name = "generic_trained_model_multi_kinase"

    frozen_folds_path = PRED_OUTPUTS_DIR / "frozen_folds_multi_kinase_sites.csv.gz"

    n_jobs_outer = 2
    workers_per_fold = 1

    params = Node2VecParams(
        dimensions=32,
        walk_length=10,
        num_walks=25,
        workers=workers_per_fold,
        p=1.0,
        q=1.0,
        window=5,
        min_count=1,
        batch_words=4,
        seed=42,
        directed=True,
    )

    print("=== Loading generic graph ===")
    graph = load_graph(GENERIC_NODES, GENERIC_EDGES, label="generic")

    print(f"Graph nodes: {len(graph.nodes):,}")
    print(f"Graph edges: {len(graph.edges):,}\n")

    print("=== Loading evaluation inputs ===")
    candidate_kinases = read_csv_auto(CANDIDATE_KINASES_OUT)
    positive_edges = read_csv_auto(frozen_folds_path)

    print(f"Candidate kinases: {len(candidate_kinases):,}")
    print(f"Frozen evaluation edges: {len(positive_edges):,}\n")

    allowed_relations = set(GENERIC_BASE_RELATIONS)

    print("=== Running leave-one-out on multi-kinase sites: generic graph ===")

    results = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=positive_edges,
        candidate_kinases=candidate_kinases,
        node2vec_params=params,
        allowed_relations=allowed_relations,
        max_folds=None,
        random_state=42,
        verbose_every=25,
        n_jobs_outer=n_jobs_outer,
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
    print("✅ Done.")


if __name__ == "__main__":
    main()