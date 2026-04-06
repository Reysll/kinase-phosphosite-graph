from __future__ import annotations

from src_prediction.config import (
    GENERIC_EDGES,
    GENERIC_NODES,
    LIVER_EDGES,
    LIVER_NODES,
    CANDIDATE_KINASES_OUT,
    PRED_OUTPUTS_DIR,
)
from src_prediction.experiment_utils import write_experiment_outputs
from src_prediction.graph_loader import load_graph
from src_prediction.io_utils import read_csv_auto
from src_prediction.leave_one_out import Node2VecParams, run_leave_one_out
from src_prediction.metrics import summarize_results, summarize_results_text
from src_prediction.relation_filters import (
    GENERIC_BASE_RELATIONS,
    PROTEIN_CORR_RELATIONS,
    SITE_CORR_RELATIONS,
)


def main() -> None:
    #Exchange the following two blocks to run the different experiments
    #First block: generic graph, no correlation edges
    #Second block: liver graph, with site correlation edges
    #To run the generic graph experiment, set:
    """
    experiment_name = "generic_cosine_parallel_baseline"
    graph_choice = "generic"
    include_site_corr = False
    include_protein_corr = False
    """

    #To run the liver graph experiment, set:

    experiment_name = "liver_cosine_with_site_corr"
    graph_choice = "liver"
    include_site_corr = True
    include_protein_corr = False


    frozen_folds_path = PRED_OUTPUTS_DIR / "frozen_folds_50.csv.gz"

    n_jobs_outer = 2
    workers_per_fold = 1

    params = Node2VecParams(
        dimensions=32,
        walk_length=10,
        num_walks=50,
        workers=workers_per_fold,
        p=1.0,
        q=1.0,
        window=5,
        min_count=1,
        batch_words=4,
        seed=42,
        directed=True,
    )

    if graph_choice == "generic":
        print("=== Loading generic graph ===")
        graph = load_graph(GENERIC_NODES, GENERIC_EDGES, label="generic")
    elif graph_choice == "liver":
        print("=== Loading liver graph ===")
        graph = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
    else:
        raise ValueError(f"Unsupported graph_choice: {graph_choice}")

    print(f"Graph nodes: {len(graph.nodes):,}")
    print(f"Graph edges: {len(graph.edges):,}\n")

    print("=== Loading evaluation inputs ===")
    candidate_kinases = read_csv_auto(CANDIDATE_KINASES_OUT)
    positive_edges = read_csv_auto(frozen_folds_path)

    print(f"Candidate kinases: {len(candidate_kinases):,}")
    print(f"Frozen evaluation edges: {len(positive_edges):,}\n")

    allowed_relations = set(GENERIC_BASE_RELATIONS)

    if include_site_corr:
        allowed_relations |= SITE_CORR_RELATIONS
    if include_protein_corr:
        allowed_relations |= PROTEIN_CORR_RELATIONS

    print("=== Experiment relation set ===")
    for rel in sorted(allowed_relations):
        print(f"  - {rel}")
    print()

    print("=== Parallel config ===")
    print(f"Outer parallel folds: {n_jobs_outer}")
    print(f"Node2Vec workers per fold: {workers_per_fold}")
    print()

    print("=== Running leave-one-out: node2vec + cosine similarity ===")

    results = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=positive_edges,
        candidate_kinases=candidate_kinases,
        node2vec_params=params,
        allowed_relations=allowed_relations,
        max_folds=None,
        random_state=42,
        verbose_every=2,
        n_jobs_outer=n_jobs_outer,
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