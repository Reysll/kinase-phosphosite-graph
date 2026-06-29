"""
Parametrized LOO runner for SLURM cluster.

Usage (called by run_loo.sh — do not run directly):
    python scripts/cluster/cluster_loo_runner.py \\
        --network liver_fc \\
        --embedding spectral \\
        --dimensions 64 \\
        --n-jobs 24 \\
        --workers 16 \\
        --trials all \\
        --output-dir /scratch/$USER/kinase_results

Supported networks:   generic | control | cancer | liver_fc
Supported embeddings: node2vec | spectral | concat
--trials all       -> frozen_trials_all.csv.gz        (all sites, recommended)
--trials multi     -> frozen_trials_multi_kinase.csv.gz (original subset)
--dimensions N     -> N eigenvectors / N2V dims per sub-embedding
                      (concat will be 2N total since it stacks both)
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# Make src_prediction importable when run as a plain script.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from src_prediction.core.config import (
    CANDIDATE_KINASES_OUT,
    GENERIC_EDGES,
    GENERIC_NODES,
    LIVER_EDGES,
    LIVER_NODES,
    PRED_OUTPUTS_DIR,
)
from src_prediction.core.experiment_utils import write_experiment_outputs
from src_prediction.core.graph_loader import load_graph
from src_prediction.core.io_utils import read_csv_auto
from src_prediction.core.metrics import summarize_results, summarize_results_text
from src_prediction.core.relation_filters import (
    GENERIC_BASE_RELATIONS,
    SITE_CORR_CANCER_RELATIONS,
    SITE_CORR_CTRL_RELATIONS,
    SITE_CORR_FC_RELATIONS,
)
from src_prediction.embeddings.embedding_strategy import (
    ConcatEmbeddingStrategy,
    Node2VecStrategy,
    SpectralEmbeddingStrategy,
)
from src_prediction.engine.leave_one_out import run_leave_one_out

NETWORK_CONFIGS = {
    "generic": {
        "nodes": GENERIC_NODES,
        "edges": GENERIC_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS),
        "label": "generic",
    },
    "control": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_CTRL_RELATIONS,
        "label": "liver",
    },
    "cancer": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_CANCER_RELATIONS,
        "label": "liver",
    },
    "liver_fc": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS,
        "label": "liver",
    },
}


def _make_strategy(
    embedding: str, dimensions: int, workers: int
) -> Node2VecStrategy | SpectralEmbeddingStrategy | ConcatEmbeddingStrategy:
    if embedding == "node2vec":
        return Node2VecStrategy(
            dimensions=dimensions,
            walk_length=10,
            num_walks=25,
            workers=workers,
            p=1.0,
            q=1.0,
            window=5,
            min_count=1,
            batch_words=10000,
            seed=42,
            directed=True,
        )
    if embedding == "spectral":
        return SpectralEmbeddingStrategy(dimensions=dimensions, directed=False, random_state=42)
    if embedding == "concat":
        return ConcatEmbeddingStrategy(
            node2vec=Node2VecStrategy(
                dimensions=dimensions,
                walk_length=10,
                num_walks=25,
                workers=workers,
                p=1.0,
                q=1.0,
                window=5,
                min_count=1,
                batch_words=10000,
                seed=42,
                directed=True,
            ),
            spectral=SpectralEmbeddingStrategy(
                dimensions=dimensions, directed=False, random_state=42
            ),
            directed=True,
        )
    raise ValueError(f"Unknown embedding: {embedding!r}")


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--network", required=True, choices=list(NETWORK_CONFIGS))
    p.add_argument("--embedding", required=True, choices=["node2vec", "spectral", "concat"])
    p.add_argument("--dimensions", type=int, default=64)
    p.add_argument("--n-jobs", type=int, default=24)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--trials", choices=["all", "multi"], default="all")
    p.add_argument("--output-dir", type=str, default="cluster_results")
    args = p.parse_args()

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    trial_tag = "allsites" if args.trials == "all" else "multikinase"
    experiment_name = f"{args.network}_{args.embedding}_dim{args.dimensions}_{trial_tag}"

    frozen_trials_path = (
        PRED_OUTPUTS_DIR / "frozen_trials_all.csv.gz"
        if args.trials == "all"
        else PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz"
    )

    cfg = NETWORK_CONFIGS[args.network]
    strategy = _make_strategy(args.embedding, args.dimensions, args.workers)

    print(f"=== Experiment: {experiment_name} ===")
    print(f"  n_jobs_outer={args.n_jobs}  workers={args.workers}  dimensions={args.dimensions}")
    print(f"  trials: {frozen_trials_path.name}")
    print(f"  output: {output_dir}")

    print(f"\n=== Loading {args.network} graph ===")
    graph = load_graph(cfg["nodes"], cfg["edges"], label=cfg["label"])
    print(f"Graph nodes: {len(graph.nodes):,}")
    print(f"Graph edges: {len(graph.edges):,}")

    candidate_kinases = read_csv_auto(CANDIDATE_KINASES_OUT)
    positive_edges = read_csv_auto(frozen_trials_path)
    print(f"\nCandidate kinases: {len(candidate_kinases):,}")
    print(f"LOO trials:        {len(positive_edges):,}")

    print(f"\n=== Running LOO: {args.network} x {args.embedding} ===")
    results = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=positive_edges,
        candidate_kinases=candidate_kinases,
        embedding_strategy=strategy,
        allowed_relations=cfg["allowed_relations"],
        max_trials=None,
        random_state=42,
        verbose_every=500,
        n_jobs_outer=args.n_jobs,
        max_negatives_per_site=50,
    )

    metrics_df = summarize_results(results)
    summary_text = summarize_results_text(results)

    paths = write_experiment_outputs(
        experiment_name=experiment_name,
        base_dir=output_dir,
        results_df=results,
        metrics_df=metrics_df,
        summary_text=summary_text,
    )

    print(summary_text)
    print(f"\nWrote: {paths['results']}")
    print(f"Wrote: {paths['metrics']}")
    print(f"Wrote: {paths['summary']}")
    print("Done.")


if __name__ == "__main__":
    main()
