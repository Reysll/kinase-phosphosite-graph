"""
Fast screening harness for model configurations.

WHY
---
Every experiment so far has cost 1 to 6 hours per configuration, which is why
only four of them have been run. The goal has changed from documenting the
current model to finding one that works, so the bottleneck is now the number of
configurations we can reject per evening. This module runs any configuration on
a fixed subsample and scores it against five pass/fail criteria.

WHY THE SUBSAMPLE IS EXACT, NOT APPROXIMATE
-------------------------------------------
`run_leave_one_out` uses `positive_edges` for three separate things: the
training matrix, the site->true-kinase map, and the trial list. Shrinking it
would shrink the TRAINING SET too, which changes the model and makes the screen
measure the wrong thing.

Instead we monkeypatch `leave_one_out._prepare_trials` (the precedent is
`null_network.py` and `random_negatives.py`) so that `positive_edges` stays the
full 6,909 edges for training while only the trial LIST is subsampled. Because
each trial fits its own model from the same full `X_train`, and the embedding is
fit once outside the loop, every retained trial produces the SAME rank it would
in a full run. The subsample is therefore exact per trial, and the only cost of
screening is statistical power, not bias.

The same holds for the predict-all pass: `positive_edges` is passed unrestricted
(as the shipped runner does), so the global model and every masked model are
identical and `target_sites` only decides which sites get scored.

This is what makes `--mode reference` free: the current model's scorecard can be
computed by restricting the EXISTING full runs to the screen subsample, and
comparing that against the full-run verdicts is a direct test of whether the
subsample preserves the conclusions.

TRIAL SUBSAMPLE DESIGN (two phase, one run)
-------------------------------------------
The frozen set is 85% rich (5,872 rich / 836 average / 201 poor). That creates a
conflict: a proportional 1,000-trial sample leaves only ~29 poor trials, but a
balanced sample has the wrong stratum mix, and the context effect lives almost
entirely in rich kinases, so a balanced sample dilutes it. Measured: a balanced
1,051-trial draw reported the Generic-vs-Liver-FC effect at |r| = 0.072 against
the full run's 0.380.

The subsample is therefore two overlapping phases, both run in one pass:

  * `prop`: 1,000 trials in the full-set proportions, used for every OVERALL
    test, so overall effect sizes stay comparable to the full runs.
  * `boost`: extra poor and average trials, bringing the union to all 201 poor
    and 350 average, used for the PER-STRATUM tests, which need power more than
    they need representativeness.

`screen_trials.csv` carries an `in_prop` flag. One LOO run covers both phases.

SITE SUBSAMPLE DESIGN
---------------------
Two overlapping sets, scored in one predict-all pass over their union:
  * `eta` set: 300 FC-valid sites drawn at random, for criterion 4. Random,
    because eta-squared is a property of the site population.
  * `kinase` set: every substrate site of the 29 kinases that have at least 5
    known substrates among the FC-valid sites, for the kinase-centric metric.
    A random 300-site draw would leave almost no kinase with 5 substrates.

WRITES
------
Only into thesis/analysis/out/screen/. Imports project code, modifies none.

Usage
-----
    python -m thesis.analysis.screen --mode subsample
    python -m thesis.analysis.screen --mode reference
    python -m thesis.analysis.screen --mode probe --network liver_fc --limit 60
    python -m thesis.analysis.screen --mode loo   --strategy spectral_orig --network liver_fc
    python -m thesis.analysis.screen --mode infer --strategy spectral_orig --network liver_fc
    python -m thesis.analysis.screen --mode scorecard --strategy spectral_orig
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from thesis.analysis.stat_tests import KEY, RANK_COL, holm, paired_test

OUT = Path("thesis/analysis/out/screen")
PRED = Path("outputs_prediction")

# ---------------------------------------------------------------- subsample design
SUBSAMPLE_SEED = 20261027
N_PROP = 1000            # phase 1: proportional, drives every OVERALL test
MIN_PER_CAT = {"poor": 201, "average": 350}   # phase 2: per-stratum power floor
N_ETA_SITES = 300
MIN_SUBSTRATES = 5       # kinase-centric metric: kinases with >= this many substrates

CAT_COL = "held_out_kinase_category"
INCAT = "adjusted_held_out_rank_in_category"
CATS = ("rich", "average", "poor")

# ---------------------------------------------------------------- pass/fail thresholds
# PROVISIONAL. Criteria 1-3 and 5 are relative, so their thresholds are forced by
# the design. Criterion 4's two numbers are judgement calls: they are set just
# beyond what node2vec already reaches (eta^2 0.990-0.992, median SD 8.2-10.2), so
# passing means being more site-aware than any model measured so far.
ALPHA = 0.05
ETA2_MAX = 0.95
MEDIAN_SD_MIN = 10.0
# Largest |rank-biserial| produced by changing the node2vec seed alone.
# source: thesis/analysis/out/seedvar_compare.txt
SEED_NOISE_CEILING = 0.373

# Reference values from the full runs, quoted in the report for orientation only.
FULL_RUN_REFS = {
    "generic": {"adj_median": 79.0, "eta2": 0.9886, "median_sd": 4.93, "n_sites": 391},
    "liver_fc": {"adj_median": 74.0, "eta2": 0.9997, "median_sd": 1.18, "n_sites": 3228},
    "control": {"adj_median": 74.0, "eta2": 0.9997, "median_sd": 1.05, "n_sites": 3228},
}

ORIG_LOO = {
    "generic": PRED / "generic_multi_kinase_spectral" / "results.csv.gz",
    "liver_fc": PRED / "liver_multi_kinase_spectral" / "results.csv.gz",
    "control": PRED / "control_multi_kinase_spectral" / "results.csv.gz",
}
ORIG_INFER = {
    "generic": PRED / "inference_all" / "generic_spectral" / "ranks.csv.gz",
    "liver_fc": PRED / "inference_all" / "liver_fc_spectral" / "ranks.csv.gz",
    "control": PRED / "inference_all" / "control_spectral" / "ranks.csv.gz",
}
ORIG_NULL = Path("thesis/analysis/out/null/null_spectral_seed1_results.csv.gz")


# ============================================================== strategies
def _spectral_orig():
    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy
    return SpectralEmbeddingStrategy(dimensions=32, directed=False, random_state=42)


def _spectral_fixed():
    """Only available on branch fix/spectral-nullspace."""
    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy
    return SpectralEmbeddingStrategy(
        dimensions=32, directed=False, random_state=42, drop_null_space=True)


def _node2vec():
    from src_prediction.embeddings.embedding_strategy import Node2VecStrategy
    # batch_words=10000 and workers=4 mirror the four Phase 1 runners exactly;
    # the dataclass defaults (4 / 8) are not comparable to any run on disk.
    return Node2VecStrategy(dimensions=32, walk_length=10, num_walks=25,
                            workers=4, p=1.0, q=1.0, window=5, min_count=1,
                            batch_words=10000, seed=42, directed=True)


STRATEGIES = {
    "spectral_orig": _spectral_orig,
    "spectral_fixed": _spectral_fixed,
    "node2vec": _node2vec,
}


# ============================================================== shared inputs
def _inputs():
    from src_prediction.core.config import (
        CANDIDATE_KINASES_OUT, GENERIC_EDGES, PRED_OUTPUTS_DIR)
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.data_prep.kinase_categories import (
        assign_kinase_categories, compute_psp_substrate_counts)

    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    counts = compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES))
    cats = assign_kinase_categories(cands["node_id"].astype(str).tolist(), counts)
    trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")
    trials = trials.reset_index(drop=True)
    trials[KEY] = np.arange(1, len(trials) + 1)
    return cands, cats, counts, trials


def scored_candidate_set() -> list[str]:
    from thesis.analysis.popularity_baseline import scored_candidate_set as s
    return s()


def _edges_for(network: str, null_seed: int | None = None):
    """Return (edges_df, allowed_relations). null_seed randomizes the FC edges."""
    from src_prediction.core.config import (
        GENERIC_EDGES, GENERIC_NODES, LIVER_EDGES, LIVER_NODES)
    from src_prediction.core.graph_loader import load_graph
    from src_prediction.core.relation_filters import (
        GENERIC_BASE_RELATIONS, SITE_CORR_CTRL_RELATIONS, SITE_CORR_FC_RELATIONS)

    if null_seed is not None:
        if network != "liver_fc":
            raise ValueError(
                "the degree-preserving null randomizes the FC correlation edges, "
                "so it is only defined for network=liver_fc")
        from thesis.analysis.null_network import build_null_edges
        edges, swap_stats = build_null_edges(null_seed)
        return edges, set(GENERIC_BASE_RELATIONS) | set(SITE_CORR_FC_RELATIONS), swap_stats

    if network == "generic":
        g = load_graph(GENERIC_NODES, GENERIC_EDGES, label="generic")
        return g.edges, set(GENERIC_BASE_RELATIONS), None
    g = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
    if network == "liver_fc":
        return g.edges, set(GENERIC_BASE_RELATIONS) | set(SITE_CORR_FC_RELATIONS), None
    if network == "control":
        return g.edges, set(GENERIC_BASE_RELATIONS) | set(SITE_CORR_CTRL_RELATIONS), None
    raise ValueError(network)


# ============================================================== 1. fixed subsamples
def build_trial_subsample(force: bool = False) -> pd.DataFrame:
    p = OUT / "screen_trials.csv"
    if p.exists() and not force:
        return pd.read_csv(p)

    OUT.mkdir(parents=True, exist_ok=True)
    _, cats, _, trials = _inputs()
    trials["category"] = trials["kinase_node_id"].astype(str).map(cats).fillna("poor")

    rng = np.random.default_rng(SUBSAMPLE_SEED)

    # phase 1: proportional allocation, largest-remainder so the total is exact
    share = trials["category"].value_counts(normalize=True)
    raw = {c: N_PROP * float(share.get(c, 0.0)) for c in CATS}
    alloc = {c: int(np.floor(v)) for c, v in raw.items()}
    for c in sorted(CATS, key=lambda c: raw[c] - alloc[c], reverse=True):
        if sum(alloc.values()) >= N_PROP:
            break
        alloc[c] += 1

    prop_idx: set[int] = set()
    for cat in CATS:
        sub = trials.loc[trials["category"] == cat, KEY].to_numpy()
        n = min(alloc[cat], len(sub))
        prop_idx |= set(rng.choice(sub, size=n, replace=False).tolist())

    # phase 2: top the thin strata up to their power floor
    union_idx = set(prop_idx)
    for cat, floor in MIN_PER_CAT.items():
        sub = trials.loc[trials["category"] == cat, KEY].to_numpy()
        have = len(prop_idx & set(sub.tolist()))
        need = min(floor, len(sub)) - have
        if need > 0:
            pool = np.array(sorted(set(sub.tolist()) - prop_idx))
            union_idx |= set(rng.choice(pool, size=min(need, len(pool)),
                                        replace=False).tolist())

    out = (trials[trials[KEY].isin(union_idx)]
           .sort_values(KEY)[[KEY, "kinase_node_id", "site_node_id", "category"]]
           .reset_index(drop=True))
    out["in_prop"] = out[KEY].isin(prop_idx)
    out.to_csv(p, index=False)
    print("[screen] trial subsample: {:,} trials, {:,} in the proportional phase".format(
        len(out), int(out["in_prop"].sum())))
    print("[screen]   union   : {}".format(out["category"].value_counts().to_dict()))
    print("[screen]   prop    : {}".format(
        out[out["in_prop"]]["category"].value_counts().to_dict()))
    return out


def build_site_subsample(force: bool = False) -> pd.DataFrame:
    p = OUT / "screen_sites.csv"
    if p.exists() and not force:
        return pd.read_csv(p)

    OUT.mkdir(parents=True, exist_ok=True)
    from src_prediction.core.config import LIVER_POSITIVE_EDGES_OUT
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    fc_sites = sorted(load_site_fc_values(LIVER_EXCEL).keys())
    scored = set(scored_candidate_set())
    pos = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    pos = pos[pos["site_node_id"].astype(str).isin(set(fc_sites))
              & pos["kinase_node_id"].astype(str).isin(scored)]

    per_kinase = pos["kinase_node_id"].astype(str).value_counts()
    rich_kinases = set(per_kinase[per_kinase >= MIN_SUBSTRATES].index)
    kinase_sites = set(
        pos[pos["kinase_node_id"].astype(str).isin(rich_kinases)]["site_node_id"].astype(str))

    rng = np.random.default_rng(SUBSAMPLE_SEED)
    eta_sites = set(rng.choice(np.array(fc_sites), size=N_ETA_SITES, replace=False).tolist())

    union = sorted(eta_sites | kinase_sites)
    out = pd.DataFrame({
        "site_node_id": union,
        "in_eta_set": [s in eta_sites for s in union],
        "in_kinase_set": [s in kinase_sites for s in union],
    })
    out.to_csv(p, index=False)
    print("[screen] site subsample: {:,} sites ({:,} eta, {:,} kinase-centric, "
          "{} kinases with >= {} substrates)".format(
              len(out), len(eta_sites), len(kinase_sites),
              len(rich_kinases), MIN_SUBSTRATES))
    return out


def _fc_valid_sites() -> set[str]:
    """The 3,228 FC-valid liver sites: the population the full Phase 3 runs used."""
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values
    return set(str(x) for x in load_site_fc_values(LIVER_EXCEL).keys())


def true_pairs_on(sites: set[str]) -> dict[str, set[str]]:
    """kinase -> its known substrate sites, restricted to *sites*."""
    from src_prediction.core.config import LIVER_POSITIVE_EDGES_OUT
    from src_prediction.core.io_utils import read_csv_auto
    pos = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    pos = pos[pos["site_node_id"].astype(str).isin(sites)]
    return {k: set(g["site_node_id"].astype(str))
            for k, g in pos.groupby(pos["kinase_node_id"].astype(str))}


# ============================================================== 2. running a config
def _patch_prepare_trials(keep: set[int]):
    """Subsample the trial LIST only; training still sees all positive_edges."""
    from src_prediction.engine import leave_one_out as loo_mod

    original = loo_mod._prepare_trials

    def patched(positive_edges, max_trials=None, random_state=42):
        t = positive_edges.copy().reset_index(drop=True)
        t[KEY] = np.arange(1, len(t) + 1)          # same indices as every full run
        t = t[t[KEY].isin(keep)].reset_index(drop=True)
        print("[screen] trial list subsampled to {:,} of {:,}".format(
            len(t), len(positive_edges)))
        return t

    loo_mod._prepare_trials = patched
    return loo_mod, original


def run_loo(strategy: str, network: str, null_seed: int | None = None,
            limit: int | None = None) -> None:
    from src_prediction.engine.leave_one_out import run_leave_one_out

    OUT.mkdir(parents=True, exist_ok=True)
    screen = build_trial_subsample()
    keep = set(screen[KEY].astype(int))
    if limit is not None:
        keep = set(sorted(keep)[:limit])

    cands, cats, _, _ = _inputs()
    from src_prediction.core.config import PRED_OUTPUTS_DIR
    from src_prediction.core.io_utils import read_csv_auto
    full_trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")

    edges, allowed, swap_stats = _edges_for(network, null_seed)
    tag = _tag(strategy, network, null_seed, limit)
    d = OUT / tag
    d.mkdir(parents=True, exist_ok=True)

    mod, original = _patch_prepare_trials(keep)
    try:
        t0 = time.time()
        res = run_leave_one_out(
            graph_edges=edges,
            positive_edges=full_trials,          # FULL set: training stays intact
            candidate_kinases=cands,
            embedding_strategy=STRATEGIES[strategy](),
            allowed_relations=allowed,
            max_trials=None,
            random_state=42,
            verbose_every=100,
            n_jobs_outer=8,
            max_negatives_per_site=50,
            kinase_categories=cats,
        )
        elapsed = time.time() - t0
    finally:
        mod._prepare_trials = original

    res.to_csv(d / "loo_results.csv.gz", index=False, compression="gzip")
    meta = {
        "tag": tag, "mode": "loo", "strategy": strategy, "network": network,
        "null_seed": null_seed, "limit": limit,
        "n_trials": int(len(res)),
        "n_rankable": int(res[RANK_COL].notna().sum()),
        "adj_median_rank": float(res[RANK_COL].median()),
        "elapsed_seconds": round(elapsed, 1),
        "seconds_per_trial": round(elapsed / max(len(res), 1), 3),
        "subsample_seed": SUBSAMPLE_SEED,
        "swap_stats": swap_stats,
    }
    (d / "loo_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("[screen] {} -> {:,} trials, median {}, {:.1f} min".format(
        tag, len(res), res[RANK_COL].median(), elapsed / 60.0))


def run_infer(strategy: str, network: str) -> None:
    from src_prediction.core.config import CANDIDATE_KINASES_OUT, LIVER_POSITIVE_EDGES_OUT
    from src_prediction.core.io_utils import read_csv_auto, write_csv_gz
    from src_prediction.engine.inference_all import run_inference_all

    OUT.mkdir(parents=True, exist_ok=True)
    sites = build_site_subsample()
    target = sorted(sites["site_node_id"].astype(str))
    # run_inference_all.py restricts positive_edges to the FULL FC-valid site set
    # before passing it, so the global model trains on 736 positives rather than
    # all 14,596. Restricting to the same 3,228 sites (NOT to `target`) keeps
    # X_train identical to the full run, so only site coverage changes.
    # Measured: passing the unrestricted 14,596 instead changes the scores
    # (median per-kinase rank SD 1.02 -> 0.29 on the same sites).
    all_pos = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    all_pos = all_pos[all_pos["site_node_id"].astype(str).isin(_fc_valid_sites())]
    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    edges, allowed, _ = _edges_for(network)

    tag = _tag(strategy, network, None, None)
    d = OUT / tag
    d.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    out = run_inference_all(
        graph_edges=edges,
        positive_edges=all_pos,
        candidate_kinases=cands,
        target_sites=target,
        embedding_strategy=STRATEGIES[strategy](),
        allowed_relations=allowed,
        max_negatives_per_site=50,
        n_jobs_outer=8,
        verbose_every=50,
    )
    elapsed = time.time() - t0

    write_csv_gz(out, d / "infer_ranks.csv.gz")
    meta = {
        "tag": tag, "mode": "infer", "strategy": strategy, "network": network,
        "target_sites": len(target),
        "sites_covered": int(out["site_node_id"].nunique()),
        "pairs_scored": int(len(out)),
        "elapsed_seconds": round(elapsed, 1),
        "subsample_seed": SUBSAMPLE_SEED,
    }
    (d / "infer_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("[screen] {} infer -> {:,}/{:,} sites, {:.1f} min".format(
        tag, out["site_node_id"].nunique(), len(target), elapsed / 60.0))


def run_eta(strategy: str, network: str) -> None:
    """Criterion 4 only, at a fraction of the cost of a full predict-all pass.

    A masked own-edge model is needed for a pair that was a positive label in
    training. On the 300-site eta set those are 47 of 124,500 pairs (0.04%), and
    eta-squared and the median per-kinase rank SD are IDENTICAL to six decimal
    places whether they are included or dropped (measured 2026-10-27). So
    criterion 4 can be scored from the global model alone, which removes the
    hundreds of per-pair model fits that dominate the predict-all runtime.

    Criterion 6 still needs the masked scores, because there the true pairs are
    the positives, so it stays on `--mode infer`.
    """
    from src_prediction.core.config import CANDIDATE_KINASES_OUT, LIVER_POSITIVE_EDGES_OUT
    from src_prediction.core.io_utils import read_csv_auto, write_csv_gz
    from src_prediction.embeddings.embeddings import build_graph_from_edges
    from src_prediction.engine.inference_all import _fit_model, _score_site
    from src_prediction.engine.leave_one_out import (
        _build_embedding_relations, _build_scoring_cache, _build_training_matrix)
    from src_prediction.data_prep.truth_utils import build_site_to_true_kinases

    OUT.mkdir(parents=True, exist_ok=True)
    sites = build_site_subsample()
    target = sorted(sites.loc[sites["in_eta_set"], "site_node_id"].astype(str))
    pos = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    pos = pos[pos["site_node_id"].astype(str).isin(_fc_valid_sites())]
    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    cand_ids = cands["node_id"].astype(str).tolist()
    edges, allowed, _ = _edges_for(network)

    t0 = time.time()
    rels = _build_embedding_relations(edges, allowed)
    strat = STRATEGIES[strategy]()
    eg = build_graph_from_edges(edges_df=edges, allowed_relations=rels,
                               directed=strat.directed)
    emb = strat.fit(eg)
    site_true = build_site_to_true_kinases(pos)
    X, y, _ = _build_training_matrix(
        positive_edges=pos, embeddings=emb, candidate_kinase_ids=cand_ids,
        site_to_true_kinases=site_true, max_negatives_per_site=50)
    cache = _build_scoring_cache(target, cand_ids, emb)
    model = _fit_model(X, y)

    rows = [_score_site(model, *cache[s], s) for s in target if s in cache]
    out = pd.concat(rows, ignore_index=True)
    out = out.sort_values(["site_node_id", "score", "kinase_node_id"],
                          ascending=[True, False, True])
    out["rank_in_site"] = out.groupby("site_node_id").cumcount() + 1
    elapsed = time.time() - t0

    tag = _tag(strategy, network, None, None)
    d = OUT / tag
    d.mkdir(parents=True, exist_ok=True)
    write_csv_gz(out.reset_index(drop=True), d / "eta_ranks.csv.gz")
    stats_ = eta_sq_and_sd(out)
    meta = {"tag": tag, "mode": "eta", "strategy": strategy, "network": network,
            "target_sites": len(target), "elapsed_seconds": round(elapsed, 1),
            "masking": "none (see run_eta docstring)", **stats_}
    (d / "eta_meta.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("[screen] {} eta -> {:,}/{:,} sites, eta^2={:.6f}, median SD={:.4f}, "
          "{:.1f} min".format(tag, out["site_node_id"].nunique(), len(target),
                              stats_["eta_squared"],
                              stats_["median_SD_rank_per_kinase"], elapsed / 60.0))


def _tag(strategy: str, network: str, null_seed: int | None, limit: int | None) -> str:
    t = "{}__{}".format(strategy, network)
    if null_seed is not None:
        t += "__null{}".format(null_seed)
    if limit is not None:
        t += "__probe{}".format(limit)
    return t


# ============================================================== 3. baselines
def popularity_ranks(keep: set[int], cats: dict[str, str]) -> pd.DataFrame:
    from thesis.analysis.popularity_baseline import build_baseline
    return build_baseline(keep, cats, counts_from="frozen",
                          candidates=scored_candidate_set())


def degree_ranks(network: str, keep: set[int], cats: dict[str, str]) -> pd.DataFrame:
    """Criterion 2: rank candidates by their degree in the embedding graph.

    Same shape as the popularity baseline, same tie-break (key DESC, id ASC),
    same adjustment (other known true kinases for the site are dropped first).
    Degree is taken on the embedding graph, which excludes `phosphorylates`, so
    it carries no direct information about the target.
    """
    from src_prediction.core.config import PRED_OUTPUTS_DIR
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.data_prep.truth_utils import build_site_to_true_kinases
    from src_prediction.embeddings.embeddings import build_graph_from_edges
    from src_prediction.engine.leave_one_out import _build_embedding_relations

    cache = OUT / "degree_{}.csv".format(network)
    if cache.exists():
        dd = pd.read_csv(cache)
        deg = dict(zip(dd["node_id"].astype(str), dd["degree"].astype(float)))
    else:
        edges, allowed, _ = _edges_for(network)
        rels = _build_embedding_relations(edges, allowed)
        eg = build_graph_from_edges(edges_df=edges, allowed_relations=rels, directed=False)
        deg = {str(n): float(d) for n, d in eg.degree()}
        pd.DataFrame({"node_id": list(deg), "degree": list(deg.values())}).to_csv(
            cache, index=False)
        print("[screen] degree cache written for {} ({:,} nodes)".format(network, len(deg)))

    trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")
    trials = trials.reset_index(drop=True)
    trials[KEY] = np.arange(1, len(trials) + 1)
    site_true = build_site_to_true_kinases(trials)

    candidates = scored_candidate_set()
    cand = np.array(candidates)
    key = np.array([deg.get(k, 0.0) for k in candidates])
    cat_of = np.array([cats.get(k, "poor") for k in candidates])

    rows = []
    for r in trials.itertuples(index=False):
        ti = int(getattr(r, KEY))
        if ti not in keep:
            continue
        k, s = str(r.kinase_node_id), str(r.site_node_id)
        if k not in set(candidates):
            continue
        others = {x for x in site_true.get(s, []) if x != k}
        m = np.array([c not in others for c in cand])
        ids, kv = cand[m], key[m]
        pos = int(np.flatnonzero(ids[np.lexsort((ids, -kv))] == k)[0]) + 1
        mc = m & (cat_of == cats.get(k, "poor"))
        ids_c, kv_c = cand[mc], key[mc]
        pos_c = int(np.flatnonzero(ids_c[np.lexsort((ids_c, -kv_c))] == k)[0]) + 1
        rows.append({KEY: ti, "true_kinase_node_id": k, "site_node_id": s,
                     RANK_COL: pos, "rank_in_category": pos_c,
                     "category": cats.get(k, "poor")})
    return pd.DataFrame(rows)


# ============================================================== 4. metrics
def eta_sq_and_sd(ranks: pd.DataFrame) -> dict:
    """eta^2 of rank_in_site ~ kinase identity, plus per-kinase rank SD."""
    y = ranks["rank_in_site"].to_numpy(float)
    grand = y.mean()
    g = ranks.groupby("kinase_node_id", observed=True)["rank_in_site"]
    means, counts, sd = g.mean(), g.size(), g.std()
    ss_between = float((counts * (means - grand) ** 2).sum())
    ss_total = float(((y - grand) ** 2).sum())
    return {
        "n_sites": int(ranks["site_node_id"].nunique()),
        "n_kinases": int(ranks["kinase_node_id"].nunique()),
        "n_pairs": int(len(ranks)),
        "eta_squared": ss_between / ss_total if ss_total else float("nan"),
        "median_SD_rank_per_kinase": float(sd.median()),
        "pct_kinases_SD_le_1": float(100 * (sd <= 1).mean()),
    }


def kinase_centric(ranks: pd.DataFrame, truth: dict[str, set[str]],
                   min_substrates: int = MIN_SUBSTRATES) -> pd.DataFrame:
    """Criterion/metric 3: rank SITES for each kinase, not kinases for each site.

    For every kinase with at least `min_substrates` known substrates inside the
    scored site pool, record where those substrates fall among all scored sites.
    The popularity baseline is constant in this direction (a kinase's substrate
    count does not vary across sites), so it cannot score above chance here.
    """
    rows = []
    for k, g in ranks.groupby("kinase_node_id", observed=True):
        sites = g["site_node_id"].astype(str).to_numpy()
        pos_set = truth.get(k, set())
        is_pos = np.array([s in pos_set for s in sites])
        n1 = int(is_pos.sum())
        if n1 < min_substrates or n1 == len(sites):
            continue
        s = g["score"].to_numpy(float)
        site_rank = stats.rankdata(-s, method="average")      # 1 = highest score
        r = stats.rankdata(s)
        n0 = len(sites) - n1
        auc = float((r[is_pos].sum() - n1 * (n1 + 1) / 2.0) / (n1 * n0))
        rows.append({
            "kinase_node_id": k,
            "n_substrates": n1,
            "n_sites_ranked": int(len(sites)),
            "median_substrate_site_rank": float(np.median(site_rank[is_pos])),
            "median_rank_pct": float(100.0 * np.median(site_rank[is_pos]) / len(sites)),
            "auc": auc,
        })
    return pd.DataFrame(rows).sort_values("auc", ascending=False).reset_index(drop=True)


# ============================================================== 5. scorecard
def _ptest(a, b, la, lb) -> dict:
    r = paired_test(np.asarray(a, float), np.asarray(b, float), la, lb)
    diff = np.asarray(a, float) - np.asarray(b, float)
    r["B_wins"] = int((diff > 0).sum())
    r["A_wins"] = int((diff < 0).sum())
    return r


def _merge(left: pd.DataFrame, lcol: str, right: pd.DataFrame, rcol: str,
           keep: set[int]) -> pd.DataFrame:
    a = (left[left[KEY].isin(keep)][[KEY, lcol, "category"]]
         if "category" in left.columns else left[left[KEY].isin(keep)][[KEY, lcol]])
    a = a.rename(columns={lcol: "a"})
    b = right[right[KEY].isin(keep)][[KEY, rcol]].rename(columns={rcol: "b"})
    return a.merge(b, on=KEY).dropna(subset=["a", "b"])


def scorecard(label: str, network: str, loo: pd.DataFrame,
              infer: pd.DataFrame | None, null_loo: pd.DataFrame | None,
              other_network_loo: pd.DataFrame | None,
              keep: set[int], cats: dict[str, str],
              prop: set[int] | None = None) -> tuple[list[dict], dict, pd.DataFrame]:
    """Returns (test rows, criterion-4/6 summary, kinase-centric table).

    `keep` is every usable trial, used for the per-stratum tests. `prop` is the
    proportional phase, used for every OVERALL test so that overall effect sizes
    stay comparable to the full runs. prop=None means keep is already
    proportional, which is the case when scoring a full run.
    """
    rows: list[dict] = []
    prop = keep if prop is None else (prop & keep)
    loo = loo[loo[KEY].isin(keep)]
    cat_of_trial = loo.set_index(KEY)[CAT_COL].to_dict() if CAT_COL in loo.columns else {}

    def strata(m: pd.DataFrame) -> pd.Series:
        if "category" in m.columns:
            return m["category"]
        return m[KEY].map(cat_of_trial)

    def pick(m: pd.DataFrame, stratum: str) -> pd.DataFrame:
        return m[m[KEY].isin(prop)] if stratum == "ALL" else m[m["category"] == stratum]

    # ---- criterion 1: popularity baseline
    pop = popularity_ranks(keep, cats)
    m = _merge(pop, RANK_COL, loo, RANK_COL, keep)
    m["category"] = strata(m)
    for stratum in ("ALL",) + CATS:
        sub = pick(m, stratum)
        if len(sub) < 10:
            continue
        r = _ptest(sub["a"], sub["b"], "popularity", "model")
        r.update({"criterion": "1 popularity", "stratum": stratum, "n": len(sub)})
        rows.append(r)

    # ---- criterion 2: degree baseline
    degb = degree_ranks(network, keep, cats)
    m2 = _merge(degb, RANK_COL, loo, RANK_COL, keep)
    m2["category"] = strata(m2)
    for stratum in ("ALL",) + CATS:
        sub = pick(m2, stratum)
        if len(sub) < 10:
            continue
        r = _ptest(sub["a"], sub["b"], "degree", "model")
        r.update({"criterion": "2 degree", "stratum": stratum, "n": len(sub)})
        rows.append(r)

    # ---- criterion 3: degree-preserving null
    if null_loo is not None:
        m3 = _merge(null_loo, RANK_COL, loo, RANK_COL, keep)
        m3["category"] = m3[KEY].map(cat_of_trial)
        for stratum in ("ALL",) + CATS:
            sub = pick(m3, stratum)
            if len(sub) < 10:
                continue
            r = _ptest(sub["a"], sub["b"], "null", "model")
            r.update({"criterion": "3 null network", "stratum": stratum, "n": len(sub)})
            rows.append(r)

    # ---- criterion 5: context effect vs the seed-noise ceiling
    if other_network_loo is not None:
        m5 = _merge(other_network_loo, RANK_COL, loo, RANK_COL, prop)
        r = _ptest(m5["a"], m5["b"], "other network", "model")
        r.update({"criterion": "5 context", "stratum": "ALL", "n": len(m5)})
        rows.append(r)

    for r in rows:
        r["label"] = label
        r["network"] = network

    # ---- criterion 4 + metric 6 from the predict-all pass
    summary = {"label": label, "network": network}
    kc = pd.DataFrame()
    if infer is not None and len(infer):
        sites = build_site_subsample()
        eta_sites = set(sites.loc[sites["in_eta_set"], "site_node_id"].astype(str))
        kin_sites = set(sites.loc[sites["in_kinase_set"], "site_node_id"].astype(str))
        e = infer[infer["site_node_id"].astype(str).isin(eta_sites)]
        if len(e):
            summary.update(eta_sq_and_sd(e))
            summary["eta_set_coverage_pct"] = round(
                100.0 * e["site_node_id"].nunique() / max(len(eta_sites), 1), 1)
        k = infer[infer["site_node_id"].astype(str).isin(kin_sites)]
        if len(k):
            kc = kinase_centric(k, true_pairs_on(kin_sites))
            if len(kc):
                summary.update({
                    "kc_n_kinases": int(len(kc)),
                    "kc_median_auc": float(kc["auc"].median()),
                    "kc_mean_auc": float(kc["auc"].mean()),
                    "kc_pct_auc_gt_half": float(100 * (kc["auc"] > 0.5).mean()),
                    "kc_median_substrate_rank_pct": float(kc["median_rank_pct"].median()),
                })
                w = stats.wilcoxon(kc["auc"] - 0.5, alternative="greater")
                summary["kc_auc_vs_chance_p"] = float(w.pvalue)
    return rows, summary, kc


def verdicts(rows: list[dict], summary: dict) -> list[dict]:
    """Apply the thresholds. One Holm family per scorecard call."""
    ps = holm([r["p_value"] for r in rows]) if rows else []
    for r, p in zip(rows, ps):
        r["p_holm"] = p

    out = []
    by = {}
    for r in rows:
        by.setdefault((r["criterion"], r["stratum"]), r)

    def get(crit, stratum="ALL"):
        return by.get((crit, stratum))

    # criterion 1 needs ALL, poor and average
    parts, detail = [], []
    for stratum in ("ALL", "poor", "average"):
        r = get("1 popularity", stratum)
        if r is None:
            parts.append(None)
            detail.append("{}: n/a".format(stratum))
            continue
        ok = r["rank_biserial"] > 0 and r["p_holm"] < ALPHA
        parts.append(ok)
        detail.append("{}: r={:+.3f} p_holm={:.2e} {}".format(
            stratum, r["rank_biserial"], r["p_holm"], "PASS" if ok else "fail"))
    out.append({"criterion": "1 beats popularity (overall + poor + average)",
                "verdict": _verdict(parts), "detail": "; ".join(detail)})

    r = get("2 degree")
    out.append({
        "criterion": "2 beats degree baseline",
        "verdict": _verdict([None if r is None else
                             (r["rank_biserial"] > 0 and r["p_holm"] < ALPHA)]),
        "detail": "n/a" if r is None else "r={:+.3f} p_holm={:.2e} median {} vs {}".format(
            r["rank_biserial"], r["p_holm"], r["median_A"], r["median_B"])})

    r = get("3 null network")
    out.append({
        "criterion": "3 beats its own degree-preserving null",
        "verdict": _verdict([None if r is None else
                             (r["rank_biserial"] > 0 and r["p_holm"] < ALPHA)]),
        "detail": "n/a (no correlation edges to randomize)" if r is None else
                  "r={:+.3f} p_holm={:.2e} median {} vs {}".format(
                      r["rank_biserial"], r["p_holm"], r["median_A"], r["median_B"])})

    if "eta_squared" in summary:
        ok = (summary["eta_squared"] < ETA2_MAX
              and summary["median_SD_rank_per_kinase"] > MEDIAN_SD_MIN)
        d = "eta^2={:.4f} (<{} needed), median per-kinase SD={:.2f} (>{} needed), {} sites".format(
            summary["eta_squared"], ETA2_MAX, summary["median_SD_rank_per_kinase"],
            MEDIAN_SD_MIN, summary["n_sites"])
    else:
        ok, d = None, "n/a (no predict-all pass)"
    out.append({"criterion": "4 site dependence", "verdict": _verdict([ok]), "detail": d})

    r = get("5 context")
    if r is None:
        ok, d = None, "n/a (needs both networks)"
    else:
        ok = abs(r["rank_biserial"]) > SEED_NOISE_CEILING and r["p_holm"] < ALPHA
        d = "|r|={:.3f} vs seed-noise ceiling {:.3f}, p_holm={:.2e}, median {} vs {}".format(
            abs(r["rank_biserial"]), SEED_NOISE_CEILING, r["p_holm"],
            r["median_A"], r["median_B"])
    out.append({"criterion": "5 context effect exceeds seed noise",
                "verdict": _verdict([ok]), "detail": d})

    if "kc_median_auc" in summary:
        ok = (summary["kc_median_auc"] > 0.5
              and summary.get("kc_auc_vs_chance_p", 1.0) < ALPHA)
        d = ("median per-kinase AUC={:.3f}, mean={:.3f}, {:.0f}% above 0.5, "
             "p={:.2e}, {} kinases".format(
                 summary["kc_median_auc"], summary["kc_mean_auc"],
                 summary["kc_pct_auc_gt_half"], summary["kc_auc_vs_chance_p"],
                 summary["kc_n_kinases"]))
    else:
        ok, d = None, "n/a (no predict-all pass)"
    out.append({"criterion": "6 kinase-centric AUC above chance (added metric)",
                "verdict": _verdict([ok]), "detail": d})
    return out


def _verdict(parts: list) -> str:
    if not parts or all(p is None for p in parts):
        return "N/A"
    if any(p is None for p in parts):
        return "PARTIAL"
    return "PASS" if all(parts) else "FAIL"


# ============================================================== 6. reference + validation
def _load(path: Path, cols: list[str]) -> pd.DataFrame | None:
    if not path.exists():
        return None
    have = pd.read_csv(path, nrows=0).columns
    return pd.read_csv(path, usecols=[c for c in cols if c in have])


def reference() -> None:
    """Scorecard the CURRENT model from the full runs already on disk, and check
    that restricting them to the screen subsample preserves every verdict."""
    OUT.mkdir(parents=True, exist_ok=True)
    screen = build_trial_subsample()
    build_site_subsample()
    keep_screen = set(screen[KEY].astype(int))
    prop_screen = set(screen.loc[screen["in_prop"], KEY].astype(int))
    _, cats, _, _ = _inputs()

    loo_cols = [KEY, "true_kinase_node_id", RANK_COL, INCAT, CAT_COL]
    loo = {n: _load(p, loo_cols) for n, p in ORIG_LOO.items()}
    null = _load(ORIG_NULL, loo_cols)
    infer = {n: _load(p, ["site_node_id", "kinase_node_id", "score", "rank_in_site"])
             for n, p in ORIG_INFER.items()}

    bad = set()
    for df in [d for d in loo.values() if d is not None] + ([null] if null is not None else []):
        bad |= set(df.loc[df[RANK_COL].isna(), KEY].astype(int))
    keep_full = set(loo["generic"][KEY].astype(int)) - bad
    keep_screen = keep_screen - bad

    L = ["SCREENING HARNESS -- REFERENCE SCORECARD (original spectral)",
         "=" * 78, "",
         "The current model, scored by the harness. Computed from the full runs",
         "already on disk, restricted to the fixed screen subsample. Because the",
         "subsample is exact per trial (see module docstring), these are the same",
         "numbers a screen RUN would produce, with no run needed.", "",
         "Subsample: {:,} trials total {}, of which {:,} are the proportional".format(
             len(screen), screen["category"].value_counts().to_dict(), len(prop_screen)),
         "phase {} that drives the overall tests. {:,} rankable.".format(
             screen[screen["in_prop"]]["category"].value_counts().to_dict(),
             len(keep_screen)),
         "Full set:  {:,} rankable trials.".format(len(keep_full)), ""]

    all_rows, all_sum, all_verd = [], [], []
    for scope, keep in (("subsample", keep_screen), ("full", keep_full)):
        for net in ("generic", "liver_fc"):
            if loo[net] is None:
                continue
            other = loo["liver_fc"] if net == "generic" else loo["generic"]
            rows, summ, kc = scorecard(
                label="spectral_orig/{}".format(scope), network=net, loo=loo[net],
                infer=infer.get(net), null_loo=null if net == "liver_fc" else None,
                other_network_loo=other, keep=keep, cats=cats,
                prop=prop_screen if scope == "subsample" else None)
            v = verdicts(rows, summ)
            for x in v:
                x.update({"scope": scope, "network": net})
            all_rows += rows
            all_sum.append(summ)
            all_verd += v
            if scope == "subsample" and len(kc):
                kc.to_csv(OUT / "reference_kinase_centric_{}.csv".format(net), index=False)

    vd = pd.DataFrame(all_verd)
    L.append("SCORECARD")
    L.append("-" * 78)
    for net in ("generic", "liver_fc"):
        for scope in ("subsample", "full"):
            s = vd[(vd["network"] == net) & (vd["scope"] == scope)]
            if s.empty:
                continue
            L.append("")
            L.append("  {}  [{}]".format(net.upper(), scope))
            for _, r in s.iterrows():
                L.append("    {:<7} {}".format(r["verdict"], r["criterion"]))
                L.append("            {}".format(r["detail"]))
    L.append("")

    L.append("VALIDATION: does the subsample reproduce the full-run verdicts?")
    L.append("-" * 78)
    piv = vd.pivot_table(index=["network", "criterion"], columns="scope",
                         values="verdict", aggfunc="first").reset_index()
    piv["agrees"] = piv.get("subsample") == piv.get("full")
    L.append(piv.to_string(index=False))
    n_ag = int(piv["agrees"].sum())
    L.append("")
    L.append("  {} of {} verdicts agree between the subsample and the full run.".format(
        n_ag, len(piv)))
    L.append("")

    L.append("KINASE-CENTRIC METRIC (added, item 3)")
    L.append("-" * 78)
    for s in all_sum:
        if "kc_median_auc" in s and s["label"].endswith("subsample"):
            L.append("  {} {}: {} kinases, median AUC {:.3f}, mean {:.3f}, "
                     "{:.0f}% above chance, median substrate rank {:.1f}% of pool".format(
                         s["label"], s["network"], s["kc_n_kinases"], s["kc_median_auc"],
                         s["kc_mean_auc"], s["kc_pct_auc_gt_half"],
                         s["kc_median_substrate_rank_pct"]))
    L.append("")
    L.append("  Site pool for this metric is the kinase set of the site subsample,")
    L.append("  which is enriched for PSP-annotated sites by construction. The")
    L.append("  negatives are therefore other kinases' substrates, which is harder")
    L.append("  than a random site pool.")
    L.append("")

    L.append("FULL-RUN REFERENCE VALUES (for orientation)")
    L.append("-" * 78)
    L.append("  generic  spectral: adj median 79.0, eta^2 0.9886 (391 sites), SD 4.93")
    L.append("  liver_fc spectral: adj median 74.0, eta^2 0.9997 (3,228 sites), SD 1.18")
    L.append("  seed-noise ceiling |r| = 0.373   source: out/seedvar_compare.txt")
    L.append("")

    pd.DataFrame(all_rows).to_csv(OUT / "reference_tests.csv", index=False)
    pd.DataFrame(all_sum).to_csv(OUT / "reference_summary.csv", index=False)
    vd.to_csv(OUT / "reference_verdicts.csv", index=False)
    piv.to_csv(OUT / "reference_validation.csv", index=False)
    txt = "\n".join(L)
    (OUT / "reference_scorecard.txt").write_text(txt, encoding="utf-8")
    print(txt)


def probe(network: str, limit: int, strategy: str = "spectral_orig") -> None:
    """Run the harness for real, confirm exactness, and fit the cost model.

    A single probe cannot separate one-time setup from per-trial cost, and
    dividing total elapsed time by a small trial count badly overestimates a
    larger run. The probe therefore runs at two nested trial counts and solves
    `elapsed = setup + marginal * n` for both terms.
    """
    n_small = max(20, limit // 3)
    runs = []
    for n in (n_small, limit):
        run_loo(strategy, network, None, n)
        tag = _tag(strategy, network, None, n)
        meta = json.loads((OUT / tag / "loo_meta.json").read_text())
        runs.append((meta["n_trials"], meta["elapsed_seconds"], tag))

    (n1, t1, _), (n2, t2, tag2) = runs
    marginal = (t2 - t1) / (n2 - n1) if n2 != n1 else float("nan")
    setup = t1 - marginal * n1

    got = pd.read_csv(OUT / tag2 / "loo_results.csv.gz")
    ref = pd.read_csv(ORIG_LOO[network], usecols=[KEY, RANK_COL])
    m = got[[KEY, RANK_COL]].merge(ref, on=KEY, suffixes=("_screen", "_full"))
    a, b = m[RANK_COL + "_screen"], m[RANK_COL + "_full"]
    same = int(((a == b) | (a.isna() & b.isna())).sum())

    screen = build_trial_subsample()
    n_screen = len(screen)

    def est(n: int) -> float:
        return (setup + marginal * n) / 60.0

    L = ["SCREENING HARNESS -- PROBE RUN AND RUNTIME ESTIMATE",
         "=" * 78, "",
         "Config: {} on {}, nested probes at {} and {} screen trials.".format(
             strategy, network, n1, n2), "",
         "EXACTNESS CHECK",
         "-" * 78,
         "  {} of {} probe trials reproduce the full run's rank exactly.".format(same, len(m)),
         "  This is the harness's correctness test: the trial list is subsampled",
         "  but training is not, so every rank must match bit for bit.", "",
         "COST MODEL (two-point fit, not a single ratio)",
         "-" * 78,
         "  {:>5} trials: {:7.1f} s".format(n1, t1),
         "  {:>5} trials: {:7.1f} s".format(n2, t2),
         "  one-time setup  {:7.1f} s  (graph load, embedding fit, training matrix)".format(
             setup),
         "  marginal cost   {:7.3f} s per trial".format(marginal), "",
         "EXTRAPOLATION TO A FULL SCREEN",
         "-" * 78,
         "  screen subsample is {:,} trials".format(n_screen),
         "  one LOO screen         ~{:4.0f} min".format(est(n_screen)),
         "  LOO for both networks  ~{:4.0f} min".format(2 * est(n_screen)),
         "  + null network LOO     ~{:4.0f} min (plus ~95 s to randomize the edges)".format(
             est(n_screen)),
         "  = criteria 1, 2, 3, 5  ~{:4.0f} min".format(3 * est(n_screen) + 1.6),
         "",
         "  Criteria 4 and 6 need a scoring pass, whose cost is driven by the",
         "  number of known true edges on the target sites, not by the trial count:",
         "    criterion 4, --mode eta     0.3 min  (global model only, 300 sites)",
         "    criterion 6, --mode infer    46 min  (662 masked fits, 612 sites)",
         "  measured 2026-10-27 on liver FC. Criterion 6 is the whole reason the",
         "  predict-all pass stays expensive, so run it as a second stage, only on",
         "  configurations that have already cleared criteria 1 to 5.",
         "",
         "  STAGE 1 (criteria 1 to 5)  ~{:4.0f} min".format(3 * est(n_screen) + 1.6 + 0.6),
         "  STAGE 2 (criterion 6)      ~  46 min",
         "",
         "  One full 6,909-trial spectral LOO run took 62 to 96 min. Dividing the",
         "  {}-trial probe by {} naively would predict {:.0f} min for the screen".format(
             n2, n2, (t2 / n2) * n_screen / 60.0),
         "  instead of {:.0f}, which is why the fit uses two points.".format(est(n_screen)),
         ""]
    txt = "\n".join(L)
    (OUT / "probe_{}_{}.txt".format(strategy, network)).write_text(txt, encoding="utf-8")
    print(txt)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True,
                    choices=["subsample", "reference", "probe", "loo", "infer",
                             "eta", "scorecard"])
    ap.add_argument("--strategy", default="spectral_orig", choices=sorted(STRATEGIES))
    ap.add_argument("--network", choices=["generic", "control", "liver_fc"])
    ap.add_argument("--null-seed", type=int, default=None)
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--force", action="store_true")
    a = ap.parse_args()

    if a.mode == "subsample":
        build_trial_subsample(force=a.force)
        build_site_subsample(force=a.force)
    elif a.mode == "reference":
        reference()
    elif a.mode == "probe":
        probe(a.network, a.limit or 60, a.strategy)
    elif a.mode == "loo":
        run_loo(a.strategy, a.network, a.null_seed, a.limit)
    elif a.mode == "infer":
        run_infer(a.strategy, a.network)
    elif a.mode == "eta":
        run_eta(a.strategy, a.network)
    else:
        scorecard_from_disk(a.strategy)


def scorecard_from_disk(strategy: str) -> None:
    """Scorecard a configuration already run by --mode loo / --mode infer."""
    OUT.mkdir(parents=True, exist_ok=True)
    screen = build_trial_subsample()
    keep = set(screen[KEY].astype(int))
    prop = set(screen.loc[screen["in_prop"], KEY].astype(int))
    _, cats, _, _ = _inputs()
    loo_cols = [KEY, "true_kinase_node_id", RANK_COL, INCAT, CAT_COL]

    def load(tag, name):
        p = OUT / tag / name
        return _load(p, loo_cols) if name.startswith("loo") else _load(
            p, ["site_node_id", "kinase_node_id", "score", "rank_in_site"])

    nets = {}
    for net in ("generic", "liver_fc", "control"):
        d = load(_tag(strategy, net, None, None), "loo_results.csv.gz")
        if d is not None:
            nets[net] = d
    if not nets:
        raise SystemExit("no screen LOO runs found for strategy {} under {}".format(
            strategy, OUT))

    bad = set()
    for d in nets.values():
        bad |= set(d.loc[d[RANK_COL].isna(), KEY].astype(int))
    keep -= bad

    rows, sums, verds, = [], [], []
    for net, d in nets.items():
        other = nets.get("liver_fc" if net == "generic" else "generic")
        nl = load(_tag(strategy, net, 1, None), "loo_results.csv.gz")
        inf = load(_tag(strategy, net, None, None), "infer_ranks.csv.gz")
        r, s, kc = scorecard(strategy, net, d, inf, nl, other, keep, cats, prop=prop)
        v = verdicts(r, s)
        for x in v:
            x.update({"strategy": strategy, "network": net})
        rows += r
        sums.append(s)
        verds += v
        if len(kc):
            kc.to_csv(OUT / "scorecard_{}_{}_kinase_centric.csv".format(strategy, net),
                      index=False)

    L = ["SCREENING SCORECARD -- {}".format(strategy), "=" * 78, "",
         "Subsample: {:,} trials, {:,} rankable, {:,} in the proportional phase.".format(
             len(screen), len(keep), len(prop & keep)), ""]
    for net in nets:
        L.append("  {}".format(net.upper()))
        for v in [x for x in verds if x["network"] == net]:
            L.append("    {:<7} {}".format(v["verdict"], v["criterion"]))
            L.append("            {}".format(v["detail"]))
        L.append("")
    pd.DataFrame(rows).to_csv(OUT / "scorecard_{}_tests.csv".format(strategy), index=False)
    pd.DataFrame(verds).to_csv(OUT / "scorecard_{}_verdicts.csv".format(strategy), index=False)
    pd.DataFrame(sums).to_csv(OUT / "scorecard_{}_summary.csv".format(strategy), index=False)
    txt = "\n".join(L)
    (OUT / "scorecard_{}.txt".format(strategy)).write_text(txt, encoding="utf-8")
    print(txt)


if __name__ == "__main__":
    main()
