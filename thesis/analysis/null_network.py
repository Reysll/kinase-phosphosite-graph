"""
Tier 2 -- Degree-preserving random-edge null network.

QUESTION
--------
Liver FC beats Generic at trial level. Is that because the correlation edges
carry liver-cancer biology, or merely because ~977 K extra edges were added to
sites with a particular degree profile? A degree-preserving null separates the
two: it adds EXACTLY the same number of correlation edges, with EXACTLY the same
per-site correlation degree sequence, but randomized partners. Any remaining
advantage of real Liver FC over the null is attributable to which sites are
connected, i.e. to the measured co-variation structure.

DESIGN
------
  * Take the Liver FC correlation edges (site_corr_fc_pos + site_corr_fc_neg).
  * Randomize them with double-edge swaps, which preserve every node's degree
    exactly while destroying the pairing. Positive and negative edge sets are
    swapped independently so each relation keeps its own degree sequence.
  * Self-loops and duplicate pairs are rejected, so the result stays a simple
    graph with the same edge count.
  * Splice the randomized edges back in place of the real ones, leaving all
    base relations untouched, then run the standard LOO engine unchanged.
  * Repeat over several seeds to get a null distribution rather than one draw.

READ-ONLY with respect to the project: reads outputs/, imports src_prediction,
and writes only into thesis/analysis/out/null/.

Usage
-----
    # time one randomization + a short LOO to extrapolate
    .venv312/Scripts/python -m thesis.analysis.null_network --probe

    # full run
    .venv312/Scripts/python -m thesis.analysis.null_network \
        --embedding spectral --seeds 1 2 3

Seeds are explicit and recorded in every output file.
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("thesis/analysis/out/null")
FC_RELATIONS = ("site_corr_fc_pos", "site_corr_fc_neg")


# --------------------------------------------------------------- randomization
def degree_preserving_swap(
    u: np.ndarray, v: np.ndarray, *, seed: int, swap_factor: int = 10
) -> tuple[np.ndarray, np.ndarray, dict]:
    """Randomize an undirected simple edge list, preserving every degree exactly.

    Repeated double-edge swap:  (a,b),(c,d) -> (a,d),(c,b)
    accepted only when it creates neither a self-loop nor a duplicate edge.
    """
    rng = np.random.default_rng(seed)
    u = u.astype(np.int64).copy()
    v = v.astype(np.int64).copy()
    m = len(u)
    n_nodes = int(max(u.max(), v.max())) + 1

    def key(a: np.ndarray | int, b: np.ndarray | int):
        lo = np.minimum(a, b)
        hi = np.maximum(a, b)
        return lo * n_nodes + hi

    present = set(key(u, v).tolist())
    target = swap_factor * m
    done = attempts = 0
    batch = 50_000
    t0 = time.time()
    while done < target and attempts < target * 20:
        i = rng.integers(0, m, batch)
        j = rng.integers(0, m, batch)
        a, b, c, d = u[i], v[i], u[j], v[j]
        ok = (i != j) & (a != d) & (c != b) & (a != c) & (b != d)
        k1 = key(a, d)
        k2 = key(c, b)
        ok &= k1 != k2
        for idx in np.flatnonzero(ok):
            ii, jj = int(i[idx]), int(j[idx])
            aa, bb, cc, dd = int(u[ii]), int(v[ii]), int(u[jj]), int(v[jj])
            n1 = min(aa, dd) * n_nodes + max(aa, dd)
            n2 = min(cc, bb) * n_nodes + max(cc, bb)
            if n1 in present or n2 in present:
                continue
            o1 = min(aa, bb) * n_nodes + max(aa, bb)
            o2 = min(cc, dd) * n_nodes + max(cc, dd)
            present.discard(o1)
            present.discard(o2)
            present.add(n1)
            present.add(n2)
            u[ii], v[ii] = aa, dd
            u[jj], v[jj] = cc, bb
            done += 1
            if done >= target:
                break
        attempts += batch
    return u, v, {
        "edges": m, "swaps_done": done, "swaps_target": target,
        "seconds": round(time.time() - t0, 1),
    }


def build_null_edges(seed: int, swap_factor: int = 10) -> tuple[pd.DataFrame, dict]:
    """Return the liver edge table with FC correlation edges randomized."""
    from src_prediction.core.config import LIVER_EDGES
    from src_prediction.core.io_utils import read_csv_auto

    liver = read_csv_auto(LIVER_EDGES)
    is_fc = liver["relation"].isin(FC_RELATIONS)
    other = liver.loc[~is_fc]
    stats_all = {}
    pieces = [other]

    for rel in FC_RELATIONS:
        sub = liver.loc[liver["relation"] == rel]
        nodes = pd.unique(pd.concat([sub["source"], sub["target"]]).astype(str))
        idx = {n: i for i, n in enumerate(nodes)}
        u = sub["source"].astype(str).map(idx).to_numpy()
        v = sub["target"].astype(str).map(idx).to_numpy()

        deg_before = np.bincount(np.concatenate([u, v]), minlength=len(nodes))
        nu, nv, st = degree_preserving_swap(u, v, seed=seed, swap_factor=swap_factor)
        deg_after = np.bincount(np.concatenate([nu, nv]), minlength=len(nodes))
        st["degree_preserved"] = bool(np.array_equal(deg_before, deg_after))
        st["frac_edges_changed"] = float(
            np.mean((nu != u) | (nv != v))
        )
        stats_all[rel] = st

        pieces.append(pd.DataFrame({
            "source": nodes[nu], "target": nodes[nv],
            "relation": rel,
            "confidence_score": np.nan,
        }))

    null_edges = pd.concat(pieces, ignore_index=True)
    return null_edges, stats_all


# --------------------------------------------------------------- LOO
def run_null_loo(embedding: str, seed: int, max_trials: int | None) -> dict:
    from src_prediction.core.config import (
        CANDIDATE_KINASES_OUT, GENERIC_EDGES, PRED_OUTPUTS_DIR,
    )
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.core.metrics import summarize_results_text
    from src_prediction.core.relation_filters import (
        GENERIC_BASE_RELATIONS, SITE_CORR_FC_RELATIONS,
    )
    from src_prediction.data_prep.kinase_categories import (
        assign_kinase_categories, compute_psp_substrate_counts,
    )
    from src_prediction.embeddings.embedding_strategy import (
        Node2VecStrategy, SpectralEmbeddingStrategy,
    )
    from src_prediction.engine.leave_one_out import run_leave_one_out

    OUT.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    null_edges, swap_stats = build_null_edges(seed)
    t_build = time.time() - t0
    print(f"[null] edges built in {t_build:.1f}s  {json.dumps(swap_stats)}", flush=True)

    strategy = (
        SpectralEmbeddingStrategy(dimensions=32, directed=False, random_state=42)
        if embedding == "spectral"
        # Settings mirror run_multi_kinase_liver.py EXACTLY. batch_words=10000 is
        # essential: the dataclass default is 4, which both slows Word2Vec by orders
        # of magnitude and changes worker-batch interleaving, so omitting it would
        # make the null incomparable to the real Liver FC node2vec run.
        else Node2VecStrategy(dimensions=32, walk_length=10, num_walks=25,
                              workers=4, p=1.0, q=1.0, window=5, min_count=1,
                              batch_words=10000, seed=42, directed=True)
    )

    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")
    cats = assign_kinase_categories(
        cands["node_id"].astype(str).tolist(),
        compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES)),
    )

    t1 = time.time()
    res = run_leave_one_out(
        graph_edges=null_edges,
        positive_edges=trials,
        candidate_kinases=cands,
        embedding_strategy=strategy,
        allowed_relations=set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS,
        max_trials=max_trials,
        random_state=42,
        verbose_every=500,
        n_jobs_outer=8,
        max_negatives_per_site=50,
        kinase_categories=cats,
    )
    t_loo = time.time() - t1

    tag = f"null_{embedding}_seed{seed}" + ("" if max_trials is None else f"_probe{max_trials}")
    res.to_csv(OUT / f"{tag}_results.csv.gz", index=False, compression="gzip")
    meta = {
        "tag": tag, "embedding": embedding, "seed": seed,
        "max_trials": max_trials, "n_trials": int(len(res)),
        "build_seconds": round(t_build, 1), "loo_seconds": round(t_loo, 1),
        "swap_stats": swap_stats,
        "adj_median_rank": float(res["adjusted_held_out_rank"].median()),
        "adj_mean_rank": float(res["adjusted_held_out_rank"].mean()),
        "n_unrankable": int(res["adjusted_held_out_rank"].isna().sum()),
    }
    (OUT / f"{tag}_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    (OUT / f"{tag}_summary.txt").write_text(summarize_results_text(res), encoding="utf-8")
    print(json.dumps(meta, indent=2), flush=True)
    return meta


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true",
                    help="time randomization + a 200-trial LOO, then extrapolate")
    ap.add_argument("--embedding", default="spectral", choices=["spectral", "node2vec"])
    ap.add_argument("--seeds", type=int, nargs="+", default=[1])
    ap.add_argument("--max-trials", type=int, default=None)
    args = ap.parse_args()

    if args.probe:
        meta = run_null_loo(args.embedding, seed=args.seeds[0], max_trials=200)
        per_trial = meta["loo_seconds"] / max(meta["n_trials"], 1)
        print("\n--- RUNTIME ESTIMATE ---")
        print(f"build (one-time per seed)   : {meta['build_seconds']:.0f}s")
        print(f"LOO setup+200 trials        : {meta['loo_seconds']:.0f}s")
        print(f"extrapolated 6,909 trials   : "
              f"~{(meta['build_seconds'] + per_trial * 6909) / 60:.0f} min per seed")
        return

    for s in args.seeds:
        run_null_loo(args.embedding, seed=s, max_trials=args.max_trials)


if __name__ == "__main__":
    main()
