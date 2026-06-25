from __future__ import annotations

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
from src_prediction.data_prep.kinase_categories import assign_kinase_categories, compute_psp_substrate_counts
from src_prediction.engine.leave_one_out import run_leave_one_out
from src_prediction.core.metrics import summarize_results, summarize_results_text
from src_prediction.core.relation_filters import GENERIC_BASE_RELATIONS, SITE_CORR_CTRL_RELATIONS


def main() -> None:
    experiment_name = "control_multi_kinase_spectral"

    frozen_trials_path = PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz"

    strategy = SpectralEmbeddingStrategy(dimensions=32, directed=False, random_state=42)

    print("=== Loading generic edges for substrate-count categories ===")
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
    print(f"Multi-kinase LOO trials: {len(positive_edges):,}\n")

    kinase_categories = assign_kinase_categories(
        candidate_kinases["node_id"].astype(str).tolist(), substrate_counts
    )

    # Control: base relations + healthy-sample raw abundance correlation edges.
    allowed_relations = set(GENERIC_BASE_RELATIONS) | SITE_CORR_CTRL_RELATIONS

    print("=== Running leave-one-out: control, spectral embedding + categories ===")
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
        kinase_categories=kinase_categories,
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
