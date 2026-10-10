"""
Phase 3 predict-all inference on a degree-preserving null network (read-only).

Mirrors src_prediction/runners/run_inference_all.py exactly (same target sites,
same masking set, same candidate kinases, same hyperparameters) but swaps the
liver edge table for the randomized null produced by null_network.py. This lets
the Result-1/Result-2 style analyses be repeated against the null.

Writes only into thesis/analysis/out/null_inference/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.null_inference --seed 1 --embedding spectral
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

OUT = Path("thesis/analysis/out/null_inference")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--embedding", default="spectral", choices=["spectral", "node2vec"])
    args = ap.parse_args()

    from src_prediction.core.config import (
        CANDIDATE_KINASES_OUT, LIVER_POSITIVE_EDGES_OUT,
    )
    from src_prediction.core.io_utils import read_csv_auto, write_csv_gz
    from src_prediction.core.relation_filters import (
        GENERIC_BASE_RELATIONS, SITE_CORR_FC_RELATIONS,
    )
    from src_prediction.embeddings.embedding_strategy import (
        Node2VecStrategy, SpectralEmbeddingStrategy,
    )
    from src_prediction.engine.inference_all import run_inference_all
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    from thesis.analysis.null_network import build_null_edges

    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"nullinfer_{args.embedding}_seed{args.seed}"

    print("=== Loading FC-valid target sites ===")
    fc_values = load_site_fc_values(LIVER_EXCEL)
    target_sites = sorted(set(fc_values.index.astype(str)))
    print(f"Target sites: {len(target_sites):,}")

    all_pos = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    positive_edges = all_pos.loc[
        all_pos["site_node_id"].astype(str).isin(set(target_sites))
    ].reset_index(drop=True)
    print(f"Masking set: {len(positive_edges):,}")

    cands = read_csv_auto(CANDIDATE_KINASES_OUT)

    print("=== Building null edges ===")
    t0 = time.time()
    null_edges, swap_stats = build_null_edges(args.seed)
    print(f"built in {time.time() - t0:.1f}s")

    strategy = (
        SpectralEmbeddingStrategy(dimensions=32, directed=False, random_state=42)
        if args.embedding == "spectral"
        else Node2VecStrategy(dimensions=32, directed=True, seed=42)
    )

    t1 = time.time()
    out = run_inference_all(
        graph_edges=null_edges,
        positive_edges=positive_edges,
        target_sites=target_sites,
        candidate_kinases=cands,
        embedding_strategy=strategy,
        allowed_relations=set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS,
        n_jobs_outer=8,
        max_negatives_per_site=50,
        verbose_every=200,
    )
    secs = time.time() - t1

    write_csv_gz(out, OUT / f"{tag}_ranks.csv.gz")
    meta = {
        "tag": tag, "embedding": args.embedding, "seed": args.seed,
        "target_sites": len(target_sites),
        "sites_covered": int(out["site_node_id"].nunique()),
        "pairs_scored": int(len(out)),
        "true_pairs_masked": int(out["is_true_kinase"].sum()),
        "inference_seconds": round(secs, 1),
        "swap_stats": swap_stats,
    }
    (OUT / f"{tag}_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
