"""
Validation runs for the spectral null-space fix (branch fix/spectral-nullspace).

WHY
---
`SpectralEmbeddingStrategy.fit` discarded ONE eigenvector, which clears the
Laplacian null space only on a connected graph. The embedding graphs have 34
(generic), 116 (control) and 122 (liver FC) components of size >= 2, so 12 to 17
of the 32 shipped dimensions are sqrt(degree) rescaled (R^2 = 1.0000), and the
resulting kinase ranking is degree-ordered at Spearman rho = -0.9995 (findings
27, 28, 30 in thesis/STATE.md; full diagnosis in thesis/EMBEDDING_DEFECT.md).

Everything measured so far concerns EMBEDDINGS. This module runs the models, so
the open question becomes answerable: does removing the degree channel change
the held-out ranks?

WHAT IS HELD FIXED
------------------
Identical to the original spectral runs in every respect except
`drop_null_space=True`:

  * frozen trials `frozen_trials_multi_kinase.csv.gz` (6,909 trials)
  * 415 scored candidate kinases, same candidate file
  * spectral d=32, directed=False, random_state=42
  * the ORIGINAL negative sampler (the alphabetical head slice), so this is not
    confounded with the random-negatives experiment
  * 50 negatives per site, n_jobs_outer=8, liblinear with balanced weights
  * the same kinase categories (poor / average / rich)

So the comparison against `outputs_prediction/*_multi_kinase_spectral/` and
`outputs_prediction/inference_all/*_spectral/` is paired on trial_index and
differs by one flag.

WRITES
------
Only into thesis/analysis/out/spectral_fixed/. `run_leave_one_out` and
`run_inference_all` return DataFrames and write nothing themselves (checked),
so no existing results file is touched.

Usage
-----
    python -m thesis.analysis.spectral_fixed_runs --mode estimate
    python -m thesis.analysis.spectral_fixed_runs --mode loo   --network generic
    python -m thesis.analysis.spectral_fixed_runs --mode loo   --network liver_fc
    python -m thesis.analysis.spectral_fixed_runs --mode infer --network control
    python -m thesis.analysis.spectral_fixed_runs --mode infer --network liver_fc
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import pandas as pd

OUT = Path("thesis/analysis/out/spectral_fixed")
DIMS = 32
EMB_SEED = 42


def _strategy():
    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy

    return SpectralEmbeddingStrategy(
        dimensions=DIMS, directed=False, random_state=EMB_SEED, drop_null_space=True)


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


def _common_inputs():
    from src_prediction.core.config import (
        CANDIDATE_KINASES_OUT, GENERIC_EDGES, PRED_OUTPUTS_DIR)
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.data_prep.kinase_categories import (
        assign_kinase_categories, compute_psp_substrate_counts)

    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    counts = compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES))
    cats = assign_kinase_categories(cands["node_id"].astype(str).tolist(), counts)
    trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")
    return cands, cats, trials


# ---------------------------------------------------------------------- LOO
def run_loo(network: str) -> None:
    from src_prediction.core.metrics import summarize_results, summarize_results_text
    from src_prediction.engine.leave_one_out import run_leave_one_out

    OUT.mkdir(parents=True, exist_ok=True)
    tag = "fixed_loo_{}_spectral".format(network)

    cands, cats, trials = _common_inputs()
    graph, allowed = _graph_for(network)
    print("[fixed] {} graph: {:,} nodes, {:,} edges".format(
        network, len(graph.nodes), len(graph.edges)))
    print("[fixed] trials: {:,}  candidates: {:,}".format(len(trials), len(cands)))

    t0 = time.time()
    results = run_leave_one_out(
        graph_edges=graph.edges,
        positive_edges=trials,
        candidate_kinases=cands,
        embedding_strategy=_strategy(),
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
    meta = {
        "tag": tag, "mode": "loo", "network": network, "embedding": "spectral",
        "drop_null_space": True, "dimensions": DIMS, "embedding_random_state": EMB_SEED,
        "negative_sampler": "original (head slice)", "max_negatives_per_site": 50,
        "n_trials": int(len(results)),
        "n_rankable": int(results["adjusted_held_out_rank"].notna().sum()),
        "elapsed_minutes": round(elapsed / 60.0, 1),
    }
    (OUT / "{}_meta.json".format(tag)).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print(txt)
    print("[fixed] done in {:.1f} min -> {}".format(elapsed / 60.0, tag))


# ------------------------------------------------------------------- Phase 3
def run_infer(network: str) -> None:
    from src_prediction.core.config import CANDIDATE_KINASES_OUT, LIVER_POSITIVE_EDGES_OUT
    from src_prediction.core.io_utils import read_csv_auto, write_csv_gz
    from src_prediction.engine.inference_all import run_inference_all
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    OUT.mkdir(parents=True, exist_ok=True)
    tag = "fixed_infer_{}_spectral".format(network)

    target_sites = sorted(load_site_fc_values(LIVER_EXCEL).keys())
    all_pos = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    graph, allowed = _graph_for(network)
    print("[fixed] target sites: {:,}".format(len(target_sites)))

    t0 = time.time()
    out = run_inference_all(
        graph_edges=graph.edges,
        positive_edges=all_pos,
        candidate_kinases=cands,
        target_sites=target_sites,
        embedding_strategy=_strategy(),
        allowed_relations=allowed,
        max_negatives_per_site=50,
        n_jobs_outer=8,
    )
    elapsed = time.time() - t0

    write_csv_gz(out, OUT / "{}_ranks.csv.gz".format(tag))
    meta = {
        "tag": tag, "mode": "infer", "network": network, "embedding": "spectral",
        "drop_null_space": True, "dimensions": DIMS, "embedding_random_state": EMB_SEED,
        "negative_sampler": "original (head slice)",
        "target_sites": len(target_sites),
        "sites_covered": int(out["site_node_id"].nunique()),
        "pairs_scored": int(len(out)),
        "elapsed_minutes": round(elapsed / 60.0, 1),
    }
    (OUT / "{}_meta.json".format(tag)).write_text(json.dumps(meta, indent=2), encoding="utf-8")
    print("[fixed] done in {:.1f} min -> {}".format(elapsed / 60.0, tag))


# ------------------------------------------------------------------ estimate
def estimate() -> None:
    """Time the embedding at the larger k, then extrapolate the four runs."""
    import networkx as nx

    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy
    from src_prediction.embeddings.embeddings import build_graph_from_edges
    from src_prediction.engine.leave_one_out import _build_embedding_relations

    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for network in ("generic", "control", "liver_fc"):
        graph, allowed = _graph_for(network)
        rels = _build_embedding_relations(graph.edges, allowed)
        eg = build_graph_from_edges(
            edges_df=graph.edges, allowed_relations=rels, directed=False)
        n_comp = nx.number_connected_components(eg)
        row = {"network": network, "nodes": eg.number_of_nodes(),
               "edges": eg.number_of_edges(), "components": n_comp}
        for drop in (False, True):
            t0 = time.time()
            SpectralEmbeddingStrategy(
                dimensions=DIMS, directed=False, random_state=EMB_SEED,
                drop_null_space=drop).fit(eg)
            row["minutes_drop_{}".format(drop)] = round((time.time() - t0) / 60.0, 2)
        row["k_requested_fixed"] = DIMS + n_comp + 16
        rows.append(row)
        print("[estimate] {}: {:,} nodes, {} components, "
              "shipped {:.2f} min, fixed {:.2f} min".format(
                  network, row["nodes"], n_comp,
                  row["minutes_drop_False"], row["minutes_drop_True"]))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "runtime_estimate.csv", index=False)

    overhead = {r["network"]: r["minutes_drop_True"] - r["minutes_drop_False"] for r in rows}
    lines = ["RUNTIME ESTIMATE: spectral null-space fix validation runs",
             "=" * 70, "",
             "Measured embedding cost (one fit per run, outside the trial loop):",
             df.to_string(index=False), "", "Extrapolation:", ""]
    lines.append("  Spectral LOO baseline (CLAUDE.md, measured): ~62 min per run")
    lines.append("  Spectral Phase 3 baseline: ~2-5 min per network")
    lines.append("")
    lines.append("  LOO generic   : ~62 min + {:.1f} min embedding overhead".format(
        overhead["generic"]))
    lines.append("  LOO liver_fc  : ~62 min + {:.1f} min embedding overhead".format(
        overhead["liver_fc"]))
    lines.append("  infer control : ~5 min  + {:.1f} min embedding overhead".format(
        overhead["control"]))
    lines.append("  infer liver_fc: ~5 min  + {:.1f} min embedding overhead".format(
        overhead["liver_fc"]))
    lines.append("")
    total = 62 * 2 + 5 * 2 + overhead["generic"] + 2 * overhead["liver_fc"] + overhead["control"]
    lines.append("  TOTAL for all four runs: ~{:.0f} min (~{:.1f} h)".format(total, total / 60.0))
    lines.append("")
    lines.append("The embedding is fit ONCE per run, before the trial loop, so the")
    lines.append("larger k adds a fixed cost per run rather than a per-trial cost.")
    txt = "\n".join(lines)
    (OUT / "runtime_estimate.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", required=True, choices=["loo", "infer", "estimate"])
    ap.add_argument("--network", choices=["generic", "control", "liver_fc"])
    a = ap.parse_args()
    if a.mode == "estimate":
        estimate()
    elif a.mode == "loo":
        run_loo(a.network)
    else:
        run_infer(a.network)


if __name__ == "__main__":
    main()
