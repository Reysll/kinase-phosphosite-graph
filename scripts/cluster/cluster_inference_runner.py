"""
Parametrized predict-all inference runner for SLURM cluster.

Usage (called by run_inference.sh — do not run directly):
    python scripts/cluster/cluster_inference_runner.py \\
        --dimensions 64 \\
        --n-jobs 24 \\
        --workers 16 \\
        --output-dir /scratch/$USER/kinase_results

Runs all 4 networks x 3 embeddings (12 combos) in sequence within a single job.
Each combo is skipped if its ranks.csv.gz already exists (safe to resume).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

import pandas as pd

from src_prediction.analysis.ksea import load_site_fc_values
from src_prediction.core.config import (
    CANDIDATE_KINASES_OUT,
    GENERIC_EDGES,
    GENERIC_NODES,
    LIVER_EDGES,
    LIVER_NODES,
    LIVER_POSITIVE_EDGES_OUT,
    PRED_OUTPUTS_DIR,
)
from src_prediction.core.graph_loader import load_graph
from src_prediction.core.io_utils import read_csv_auto, write_csv_gz, write_text
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
from src_prediction.engine.inference_all import run_inference_all

LIVER_EXCEL = Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")

NETWORK_CONFIGS = {
    "generic": {
        "nodes": GENERIC_NODES,
        "edges": GENERIC_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS),
    },
    "control": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_CTRL_RELATIONS,
    },
    "cancer": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_CANCER_RELATIONS,
    },
    "liver_fc": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS,
    },
}


def _make_strategies(dimensions: int, workers: int) -> dict:
    return {
        "node2vec": Node2VecStrategy(
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
        "spectral": SpectralEmbeddingStrategy(
            dimensions=dimensions, directed=False, random_state=42
        ),
        "concat": ConcatEmbeddingStrategy(
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
        ),
    }


def _write_summary(out: pd.DataFrame, target_sites: list, network: str, embedding: str) -> str:
    n_sites_total = len(target_sites)
    n_sites_covered = out["site_node_id"].nunique()
    n_true_edges = int(out["is_true_kinase"].sum())
    lines = [
        f"Inference-all summary -- network={network}, embedding={embedding}",
        "=" * 60,
        f"Target sites (FC-valid): {n_sites_total:,}",
        f"Sites with >=1 scored kinase: {n_sites_covered:,} "
        f"({n_sites_covered / n_sites_total:.1%})",
        f"Total (kinase, site) pairs scored: {len(out):,}",
        f"Pairs that are known true KSAs (own-edge masked): {n_true_edges:,}",
        "",
        "Per-kinase average rank (top 20 by lowest mean rank):",
    ]
    avg_rank = (
        out.groupby("kinase_node_id")["rank_in_site"]
        .agg(["mean", "median", "count"])
        .sort_values("mean")
        .head(20)
    )
    for kinase, row in avg_rank.iterrows():
        lines.append(
            f"  {kinase:<20}  mean={row['mean']:.1f}  median={row['median']:.0f}"
            f"  n={int(row['count'])}"
        )
    return "\n".join(lines)


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--dimensions", type=int, default=64)
    p.add_argument("--n-jobs", type=int, default=24)
    p.add_argument("--workers", type=int, default=16)
    p.add_argument("--output-dir", type=str, default="cluster_results")
    args = p.parse_args()

    output_dir = Path(args.output_dir)
    inference_dir = output_dir / "inference_all"
    inference_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"=== Predict-all inference: dim={args.dimensions}  n_jobs={args.n_jobs}"
        f"  workers={args.workers} ==="
    )

    print("\n=== Loading FC-valid target sites ===")
    fc_values = load_site_fc_values(LIVER_EXCEL)
    target_sites = sorted(set(fc_values.index.astype(str)))
    print(f"Target sites: {len(target_sites):,}")

    all_positive_edges = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    target_site_set = set(target_sites)
    positive_edges = all_positive_edges.loc[
        all_positive_edges["site_node_id"].astype(str).isin(target_site_set)
    ].reset_index(drop=True)
    print(f"Restricted positive edges (masking set): {len(positive_edges):,}")

    candidate_kinases = read_csv_auto(CANDIDATE_KINASES_OUT)
    strategies = _make_strategies(args.dimensions, args.workers)
    combos = [(net, emb) for net in NETWORK_CONFIGS for emb in strategies]

    for i, (network, embedding) in enumerate(combos, 1):
        out_dir = inference_dir / f"{network}_{embedding}"
        ranks_path = out_dir / "ranks.csv.gz"
        if ranks_path.exists():
            print(f"[{i}/{len(combos)}] {network} x {embedding} -- already exists, skipping.")
            continue

        print(f"\n{'=' * 70}")
        print(f"[{i}/{len(combos)}] RUNNING: network={network}  embedding={embedding}")
        print("=" * 70)

        cfg = NETWORK_CONFIGS[network]
        graph = load_graph(cfg["nodes"], cfg["edges"], label=network)

        out = run_inference_all(
            graph_edges=graph.edges,
            positive_edges=positive_edges,
            target_sites=target_sites,
            candidate_kinases=candidate_kinases,
            embedding_strategy=strategies[embedding],
            allowed_relations=cfg["allowed_relations"],
            n_jobs_outer=args.n_jobs,
            max_negatives_per_site=50,
            verbose_every=200,
        )

        out_dir.mkdir(parents=True, exist_ok=True)
        write_csv_gz(out, ranks_path)
        summary = _write_summary(out, target_sites, network, embedding)
        write_text(summary, out_dir / "summary.txt")
        print(summary)
        print(f"\nWrote: {ranks_path}")

    print("\nAll combos complete.")


if __name__ == "__main__":
    main()
