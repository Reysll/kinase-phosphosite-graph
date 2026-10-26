"""
Random-negatives experiment (thesis-side; no project code is modified).

WHY
---
`src_prediction/data_prep/negative_sampling.py` takes `negatives[:50]`, a head
slice of an alphabetically sorted candidate list, so only 55 of 420 kinases ever
serve as negatives (finding 23). The existing ever-negative vs never-negative
comparison tests only a LOCAL effect, whether those 55 are pushed down relative
to the rest. It cannot see the GLOBAL effect: if the classifier only ever
separates positives from the same early-alphabet kinases, it fits a direction in
feature space that moves every kinase's score, which is invisible to a
between-group comparison and is a candidate mechanism for the near-fixed global
kinase ordering (eta^2 = 0.9997) and for the strength of the popularity baseline.

Testing the global effect requires rerunning with a uniformly random sampler.
Because the frozen trials, candidate set, embeddings and classifier are all held
fixed and only the negative selection changes, this is a valid paired comparison
against the original runs.

HOW
---
`sample_negative_pairs_for_site` is imported by name into
`src_prediction.engine.leave_one_out` and `src_prediction.engine.model_scoring`,
so the replacement is bound on those module objects at runtime, the same
approach `null_network.py` uses to swap the edge table. Nothing on disk under
`src_prediction/` is touched.

`src_prediction.engine.inference_all` reuses `_build_training_matrix` from
`leave_one_out`, so patching `leave_one_out` covers the Phase 3 path too. Both
engines parallelize with threads, which share module state, so the patch holds
inside workers.

The replacement draws `max_negatives` kinases uniformly without replacement from
the non-true candidates. Its RNG is seeded per site from a hash of
(sampler seed, site id), so the draw is reproducible and independent of the
order in which threads happen to reach each site.

Everything else matches the original spectral runs: same frozen trials, 415
scored candidates, liblinear with balanced weights, spectral d=32 with
random_state=42, n_jobs_outer=8, 50 negatives per site.

Writes only into thesis/analysis/out/random_negatives/.

Usage
-----
    python -m thesis.analysis.random_negatives --mode loo   --network generic  --seed 1
    python -m thesis.analysis.random_negatives --mode loo   --network liver_fc --seed 1
    python -m thesis.analysis.random_negatives --mode infer --network generic  --seed 1
    python -m thesis.analysis.random_negatives --mode estimate
"""

from __future__ import annotations

import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("thesis/analysis/out/random_negatives")
NEG_COLS = ["kinase_node_id", "site_node_id", "label"]


# --------------------------------------------------------------- the sampler
def make_random_sampler(seed: int, usage_log: dict):
    """Uniform random negatives, seeded per site so draws are order-independent."""

    def sampler(site_node_id, candidate_kinase_ids, true_kinases_for_site,
                max_negatives=None):
        negatives = [k for k in candidate_kinase_ids if k not in true_kinases_for_site]
        if max_negatives is not None and len(negatives) > max_negatives:
            key = "{}|{}".format(seed, site_node_id).encode()
            h = hashlib.blake2b(key, digest_size=8).digest()
            rng = np.random.default_rng(int.from_bytes(h, "big"))
            idx = rng.choice(len(negatives), size=max_negatives, replace=False)
            negatives = [negatives[i] for i in sorted(idx)]
        for k in negatives:
            usage_log[k] = usage_log.get(k, 0) + 1
        return pd.DataFrame({
            "kinase_node_id": negatives,
            "site_node_id": [site_node_id] * len(negatives),
            "label": [0] * len(negatives),
        })

    return sampler


def install(sampler) -> list:
    """Bind the replacement on every module that imported the original."""
    import importlib
    targets = [
        "src_prediction.engine.leave_one_out",
        "src_prediction.engine.model_scoring",
        "src_prediction.data_prep.negative_sampling",
    ]
    restores = []
    for name in targets:
        mod = importlib.import_module(name)
        if hasattr(mod, "sample_negative_pairs_for_site"):
            restores.append((mod, mod.sample_negative_pairs_for_site))
            mod.sample_negative_pairs_for_site = sampler
            print("[patch] {}.sample_negative_pairs_for_site".format(name))
    if not restores:
        raise RuntimeError("no patch target found")
    return restores


def verify_patch(sampler, seed: int) -> None:
    """Fail loudly if the replacement is not actually being used.

    The self-test uses a THROWAWAY sampler with its own usage log, so the
    synthetic kinase ids below never reach the run's usage record. An earlier
    version called `sampler` directly, which logged 88 synthetic ids and made
    the recorded distinct-kinase count exceed the 420 real candidates.
    """
    from src_prediction.engine import leave_one_out as loo
    assert loo.sample_negative_pairs_for_site is sampler, "patch not bound on leave_one_out"
    probe = make_random_sampler(seed, {})
    cands = ["PROTEIN:K{:03d}".format(i) for i in range(200)]
    a = probe("SITE:X-S1", cands, set(), max_negatives=50)
    b = probe("SITE:X-S1", cands, set(), max_negatives=50)
    c = probe("SITE:Y-S2", cands, set(), max_negatives=50)
    assert len(a) == 50, "wrong count"
    assert list(a["kinase_node_id"]) == list(b["kinase_node_id"]), "not reproducible"
    assert list(a["kinase_node_id"]) != list(c["kinase_node_id"]), "same draw for two sites"
    head = cands[:50]
    assert list(a["kinase_node_id"]) != head, "draw equals the alphabetical head slice"
    print("[patch] verified: 50 drawn, reproducible per site, differs per site, "
          "not the head slice")


# ------------------------------------------------------------------- shared
def _common_inputs():
    from src_prediction.core.config import CANDIDATE_KINASES_OUT, GENERIC_EDGES, PRED_OUTPUTS_DIR
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.data_prep.kinase_categories import (
        assign_kinase_categories, compute_psp_substrate_counts)

    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    generic_edges = read_csv_auto(GENERIC_EDGES)
    counts = compute_psp_substrate_counts(generic_edges)
    cats = assign_kinase_categories(cands["node_id"].astype(str).tolist(), counts)
    trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")
    return cands, cats, trials


def _graph_for(network: str):
    from src_prediction.core.config import (
        GENERIC_EDGES, GENERIC_NODES, LIVER_EDGES, LIVER_NODES)
    from src_prediction.core.graph_loader import load_graph
    from src_prediction.core.relation_filters import (
        GENERIC_BASE_RELATIONS, SITE_CORR_CTRL_RELATIONS, SITE_CORR_FC_RELATIONS)

    if network == "generic":
        g = load_graph(GENERIC_NODES, GENERIC_EDGES, label="generic")
        return g, set(GENERIC_BASE_RELATIONS)
    g = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
    if network == "liver_fc":
        return g, set(GENERIC_BASE_RELATIONS) | set(SITE_CORR_FC_RELATIONS)
    if network == "control":
        return g, set(GENERIC_BASE_RELATIONS) | set(SITE_CORR_CTRL_RELATIONS)
    raise ValueError(network)


# ---------------------------------------------------------------------- LOO
def run_loo(network: str, seed: int) -> None:
    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy
    from src_prediction.engine.leave_one_out import run_leave_one_out
    from src_prediction.core.metrics import summarize_results, summarize_results_text

    OUT.mkdir(parents=True, exist_ok=True)
    tag = "randneg_loo_{}_spectral_seed{}".format(network, seed)

    usage: dict = {}
    sampler = make_random_sampler(seed, usage)
    install(sampler)
    verify_patch(sampler, seed)

    cands, cats, trials = _common_inputs()
    graph, allowed = _graph_for(network)
    print("[randneg] {} graph: {:,} nodes, {:,} edges".format(
        network, len(graph.nodes), len(graph.edges)))
    print("[randneg] trials: {:,}  candidates: {:,}".format(len(trials), len(cands)))

    t0 = time.time()
    results = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=trials,
        candidate_kinases=cands,
        embedding_strategy=SpectralEmbeddingStrategy(
            dimensions=32, directed=False, random_state=42),
        allowed_relations=allowed,
        max_trials=None,
        random_state=42,
        verbose_every=250,
        n_jobs_outer=8,
        max_negatives_per_site=50,
        kinase_categories=cats,
    )
    elapsed = time.time() - t0

    results.to_csv(OUT / "{}_results.csv.gz".format(tag), index=False, compression="gzip")
    summarize_results(results).to_csv(OUT / "{}_metrics.csv".format(tag), index=False)
    txt = summarize_results_text(results)
    (OUT / "{}_summary.txt".format(tag)).write_text(txt, encoding="utf-8")

    use = pd.Series(usage).sort_values(ascending=False)
    use.to_csv(OUT / "{}_negative_usage.csv".format(tag), header=["n_uses"])
    meta = {
        "tag": tag, "mode": "loo", "network": network, "sampler_seed": seed,
        "embedding": "spectral", "dimensions": 32, "embedding_random_state": 42,
        "max_negatives_per_site": 50, "n_trials": int(len(results)),
        "n_rankable": int(results["adjusted_held_out_rank"].notna().sum()),
        "distinct_kinases_used_as_negative": int(len(use)),
        "total_negative_pairs": int(use.sum()),
        "elapsed_minutes": round(elapsed / 60.0, 1),
    }
    (OUT / "{}_meta.json".format(tag)).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(txt)
    print("[randneg] distinct kinases used as negatives: {}".format(len(use)))
    print("[randneg] done in {:.1f} min -> {}".format(elapsed / 60.0, tag))


# ------------------------------------------------------------------- Phase 3
def run_infer(network: str, seed: int) -> None:
    from src_prediction.core.config import (
        CANDIDATE_KINASES_OUT, LIVER_POSITIVE_EDGES_OUT)
    from src_prediction.core.io_utils import read_csv_auto, write_csv_gz
    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy
    from src_prediction.engine.inference_all import run_inference_all
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    OUT.mkdir(parents=True, exist_ok=True)
    tag = "randneg_infer_{}_spectral_seed{}".format(network, seed)

    usage: dict = {}
    sampler = make_random_sampler(seed, usage)
    install(sampler)
    verify_patch(sampler, seed)

    fc = load_site_fc_values(LIVER_EXCEL)
    target_sites = sorted(fc.keys())
    all_pos = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    graph, allowed = _graph_for(network)
    print("[randneg] target sites: {:,}".format(len(target_sites)))

    t0 = time.time()
    out = run_inference_all(
        graph_edges=graph.edges,
        positive_edges=all_pos,
        candidate_kinases=cands,
        target_sites=target_sites,
        embedding_strategy=SpectralEmbeddingStrategy(
            dimensions=32, directed=False, random_state=42),
        allowed_relations=allowed,
        max_negatives_per_site=50,
        n_jobs_outer=8,
    )
    elapsed = time.time() - t0

    write_csv_gz(out, OUT / "{}_ranks.csv.gz".format(tag))
    use = pd.Series(usage).sort_values(ascending=False)
    use.to_csv(OUT / "{}_negative_usage.csv".format(tag), header=["n_uses"])
    meta = {
        "tag": tag, "mode": "infer", "network": network, "sampler_seed": seed,
        "embedding": "spectral", "target_sites": len(target_sites),
        "sites_covered": int(out["site_node_id"].nunique()),
        "pairs_scored": int(len(out)),
        "distinct_kinases_used_as_negative": int(len(use)),
        "elapsed_minutes": round(elapsed / 60.0, 1),
    }
    (OUT / "{}_meta.json".format(tag)).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("[randneg] distinct kinases used as negatives: {}".format(len(use)))
    print("[randneg] done in {:.1f} min -> {}".format(elapsed / 60.0, tag))


# ------------------------------------------------------------------ estimate
def estimate() -> None:
    print("RUNTIME ESTIMATE (from measured runtimes in CLAUDE.md and run metadata)")
    print("=" * 74)
    print("  Spectral LOO, 6,909 trials      : ~62 min per run")
    print("  Phase 3 spectral inference      : ~2-5 min per network")
    print("     (nullinfer_spectral_seed1 recorded inference_seconds = 70.9")
    print("      plus embedding and training-matrix build)")
    print()
    print("  Step (a) LOO Generic + Liver FC, seed 1   : ~2.1 h")
    print("  Step (b) Phase 3 Generic/Control/Liver FC : ~0.3 h")
    print("  Step (c) LOO x2 networks x seeds 2 and 3  : ~4.2 h")
    print("  ----------------------------------------------------")
    print("  (a) + (b)              : ~2.4 h")
    print("  (a) + (b) + (c)        : ~6.6 h")
    print()
    print("  Negative selection itself is a trivial share of runtime; the")
    print("  trial loop dominates, as established for the node2vec runs.")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["loo", "infer", "estimate"])
    ap.add_argument("--network", choices=["generic", "control", "liver_fc"])
    ap.add_argument("--seed", type=int, default=1)
    a = ap.parse_args()
    if a.mode == "estimate":
        estimate()
    elif a.mode == "loo":
        run_loo(a.network, a.seed)
    else:
        run_infer(a.network, a.seed)


if __name__ == "__main__":
    main()
