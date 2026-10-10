"""
node2vec seed-variance baseline (read-only w.r.t. the project).

WHY
---
node2vec is stochastic: random walks plus Word2Vec initialisation. Every
reported Phase 1 number comes from a SINGLE seed (42). Before any Liver-FC-vs-
Generic difference can be called real, we need to know how much the metric
moves when only the seed changes. This reruns the real Generic and Liver FC
node2vec LOO at new seeds, everything else identical to the project's own
run_multi_kinase_{generic,liver}.py.

Spectral is deterministic (ARPACK with a fixed starting vector), so it needs no
seed-variance baseline.

Writes only into thesis/analysis/out/seed_variance/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.seed_variance --network generic --seed 7
    .venv312/Scripts/python -m thesis.analysis.seed_variance --network liver_fc --seed 7
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

OUT = Path("thesis/analysis/out/seed_variance")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--network", required=True, choices=["generic", "liver_fc"])
    ap.add_argument("--seed", type=int, required=True)
    args = ap.parse_args()

    from src_prediction.core.config import (
        CANDIDATE_KINASES_OUT, GENERIC_EDGES, GENERIC_NODES,
        LIVER_EDGES, LIVER_NODES, PRED_OUTPUTS_DIR,
    )
    from src_prediction.core.graph_loader import load_graph
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.core.metrics import summarize_results_text
    from src_prediction.core.relation_filters import (
        GENERIC_BASE_RELATIONS, SITE_CORR_FC_RELATIONS,
    )
    from src_prediction.data_prep.kinase_categories import (
        assign_kinase_categories, compute_psp_substrate_counts,
    )
    from src_prediction.embeddings.embedding_strategy import Node2VecStrategy
    from src_prediction.engine.leave_one_out import run_leave_one_out

    OUT.mkdir(parents=True, exist_ok=True)
    tag = f"seedvar_{args.network}_node2vec_seed{args.seed}"

    if args.network == "generic":
        graph = load_graph(GENERIC_NODES, GENERIC_EDGES, label="generic")
        allowed = set(GENERIC_BASE_RELATIONS)
    else:
        graph = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
        allowed = set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS

    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")
    cats = assign_kinase_categories(
        cands["node_id"].astype(str).tolist(),
        compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES)),
    )

    # Identical to run_multi_kinase_{generic,liver}.py in EVERY argument except
    # `seed`, which is the variable under study. batch_words=10000 must be given
    # explicitly: the dataclass default is 4, which would change Word2Vec batching
    # and make these runs incomparable to the original seed-42 results.
    strategy = Node2VecStrategy(
        dimensions=32, walk_length=10, num_walks=25, workers=4,
        p=1.0, q=1.0, window=5, min_count=1, batch_words=10000,
        seed=args.seed, directed=True,
    )

    t0 = time.time()
    res = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=trials,
        candidate_kinases=cands,
        embedding_strategy=strategy,
        allowed_relations=allowed,
        max_trials=None,
        random_state=42,
        verbose_every=500,
        n_jobs_outer=8,
        max_negatives_per_site=50,
        kinase_categories=cats,
    )
    secs = time.time() - t0

    res.to_csv(OUT / f"{tag}_results.csv.gz", index=False, compression="gzip")
    meta = {
        "tag": tag, "network": args.network, "embedding": "node2vec",
        "node2vec_seed": args.seed, "n_trials": int(len(res)),
        "loo_seconds": round(secs, 1),
        "adj_median_rank": float(res["adjusted_held_out_rank"].median()),
        "adj_mean_rank": float(res["adjusted_held_out_rank"].mean()),
        "n_unrankable": int(res["adjusted_held_out_rank"].isna().sum()),
    }
    (OUT / f"{tag}_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    (OUT / f"{tag}_summary.txt").write_text(summarize_results_text(res), encoding="utf-8")
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
