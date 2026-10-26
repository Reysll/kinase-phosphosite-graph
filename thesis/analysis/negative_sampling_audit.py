"""
Audit of negative-pair selection in the LOO training set (read-only).

CONCERN
-------
src_prediction/data_prep/negative_sampling.py builds negatives for a site as

    negatives = [k for k in candidate_kinase_ids if k not in true_kinases_for_site]
    if max_negatives is not None:
        negatives = negatives[:max_negatives]

This is a HEAD SLICE of the candidate list, not a random sample. Whether that
matters depends entirely on the order of the candidate list, which this script
establishes, and then on which kinases consequently ever appear as negatives.

The training caller passes max_negatives_per_site = 50
(src_prediction/engine/model_scoring.py).

This script reproduces the selection exactly, without importing or modifying
project code, and measures the consequences.

Reports:
  (a) the order of outputs_prediction/candidate_kinases.csv.gz
  (b) distinct kinases ever used as negatives, and the count per kinase
  (c) how many of the 415 scored kinases never appear as a negative
  (d) whether never-negative kinases rank systematically higher in the
      existing results, by LOO adjusted rank and predict-all rank

Inputs (read-only):
    outputs_prediction/candidate_kinases.csv.gz
    outputs_prediction/liver_positive_kinase_site_edges.csv.gz
    outputs_prediction/frozen_trials_multi_kinase.csv.gz
    outputs_prediction/*multi_kinase*spectral*/ or results files for LOO ranks
    outputs_prediction/inference_all/{control,liver_fc}_spectral/ranks.csv.gz

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.negative_sampling_audit
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PRED = Path("outputs_prediction")
OUT = Path("thesis/analysis/out")
INFER = PRED / "inference_all"
MAX_NEG = 50


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L = []
    add = L.append
    add("NEGATIVE-SAMPLING AUDIT")
    add("=" * 78)
    add("")
    add("Reproduces src_prediction/data_prep/negative_sampling.py exactly.")
    add("Read-only: no project code, data or result is modified.")
    add("")

    # ---------------------------------------------------------------- (a)
    cand = pd.read_csv(PRED / "candidate_kinases.csv.gz")
    order = cand["node_id"].astype(str).tolist()
    add("(a) ORDER OF THE CANDIDATE LIST")
    add("-" * 78)
    add("  file   : outputs_prediction/candidate_kinases.csv.gz")
    add("  built  : src_prediction/data_prep/evaluation_set.py, get_candidate_kinases()")
    add("  rows   : {}".format(len(order)))
    is_sorted = order == sorted(order)
    add("  sorted ascending by node_id: {}".format(is_sorted))
    add("  the builder calls .sort_values(\"node_id\"), so the order is")
    add("  ALPHABETICAL BY GENE SYMBOL, not random and not by any model quantity.")
    add("")
    add("  first 25: " + ", ".join(k.split(":")[-1] for k in order[:25]))
    add("  last 10 : " + ", ".join(k.split(":")[-1] for k in order[-10:]))
    add("")

    # --------------------------------------------- site -> true kinases
    pos = pd.read_csv(PRED / "liver_positive_kinase_site_edges.csv.gz")
    kcol = "kinase_node_id" if "kinase_node_id" in pos.columns else pos.columns[0]
    scol = "site_node_id" if "site_node_id" in pos.columns else pos.columns[1]
    site_true = pos.groupby(scol)[kcol].apply(lambda s: set(s.astype(str))).to_dict()

    trials = pd.read_csv(PRED / "frozen_trials_multi_kinase.csv.gz")
    tsite = "site_node_id" if "site_node_id" in trials.columns else scol
    train_sites = sorted(set(trials[tsite].astype(str)))

    add("TRAINING SITES CONSIDERED")
    add("-" * 78)
    add("  positive edge file rows      : {:,}".format(len(pos)))
    add("  frozen trials                : {:,}".format(len(trials)))
    add("  distinct sites in frozen set : {:,}".format(len(train_sites)))
    add("  max_negatives_per_site       : {} (model_scoring.py default)".format(MAX_NEG))
    add("")

    # ---------------------------------------------------------------- (b)
    counts = {}
    for site in train_sites:
        true_set = site_true.get(site, set())
        negs = [k for k in order if k not in true_set][:MAX_NEG]
        for k in negs:
            counts[k] = counts.get(k, 0) + 1

    cser = pd.Series(counts).sort_values(ascending=False)
    add("(b) KINASES EVER USED AS NEGATIVES")
    add("-" * 78)
    add("  distinct kinases ever a negative : {} of {} in the candidate list".format(
        len(cser), len(order)))
    add("  total negative pairs generated   : {:,}".format(int(cser.sum())))
    add("")
    add("  count per kinase (all of them, descending):")
    for k, v in cser.items():
        add("    {:<12s} {:>6,}".format(k.split(":")[-1], int(v)))
    add("")
    add("  position of these kinases in the candidate list:")
    pos_idx = sorted(order.index(k) for k in cser.index)
    add("    min index {}, max index {}, out of {}".format(
        pos_idx[0], pos_idx[-1], len(order) - 1))
    add("    they occupy the FIRST {} slots of the alphabetical list, with the".format(
        pos_idx[-1] + 1))
    add("    few extra slots beyond 50 explained by sites whose true kinases")
    add("    fall early alphabetically and are therefore skipped.")
    add("")

    # ---------------------------------------------------------------- (c)
    # scored kinases = those appearing in the predict-all output
    scored = None
    ctrl_ranks = INFER / "control_spectral" / "ranks.csv.gz"
    if ctrl_ranks.exists():
        dfr = pd.read_csv(ctrl_ranks, usecols=["kinase_node_id", "rank_in_site"])
        scored = sorted(set(dfr["kinase_node_id"].astype(str)))
    add("(c) SCORED KINASES THAT NEVER APPEAR AS A NEGATIVE")
    add("-" * 78)
    if scored is None:
        add("  predict-all output missing; cannot determine the scored set")
    else:
        ever = set(cser.index)
        never = [k for k in scored if k not in ever]
        add("  kinases scored in every run : {}".format(len(scored)))
        add("  of those, ever a negative   : {}".format(len(set(scored) & ever)))
        add("  of those, NEVER a negative  : {}  ({:.1f}%)".format(
            len(never), 100.0 * len(never) / len(scored)))
        add("")
        add("  The classifier therefore never sees an explicit negative example")
        add("  for {} of the {} kinases it is asked to rank.".format(
            len(never), len(scored)))
        add("")

        # ------------------------------------------------------------ (d)
        add("(d) DO NEVER-NEGATIVE KINASES RANK HIGHER?")
        add("-" * 78)
        rows = []

        # predict-all median rank per kinase, both context networks
        for net in ("control", "liver_fc"):
            p = INFER / "{}_spectral".format(net) / "ranks.csv.gz"
            if not p.exists():
                continue
            d = pd.read_csv(p, usecols=["kinase_node_id", "rank_in_site"])
            med = d.groupby("kinase_node_id", observed=True)["rank_in_site"].median()
            a = med[[k for k in med.index if k in ever]].to_numpy()
            b = med[[k for k in med.index if k not in ever]].to_numpy()
            if len(a) and len(b):
                u, pv = stats.mannwhitneyu(a, b, alternative="two-sided")
                # rank-biserial from U
                rb = 1.0 - 2.0 * u / (len(a) * len(b))
                rows.append({
                    "source": "predict-all {} spectral".format(net),
                    "metric": "median rank_in_site",
                    "n_ever_negative": len(a), "median_ever": float(np.median(a)),
                    "n_never_negative": len(b), "median_never": float(np.median(b)),
                    "mannwhitney_p": float(pv), "rank_biserial": float(rb),
                })

        # LOO per-kinase median adjusted rank, from the raw run outputs
        for run in ("generic_multi_kinase_spectral", "control_multi_kinase_spectral",
                    "liver_multi_kinase_spectral"):
            rp = PRED / run / "results.csv.gz"
            if not rp.exists():
                continue
            r = pd.read_csv(rp, usecols=["true_kinase_node_id",
                                         "adjusted_held_out_rank",
                                         "adjusted_held_out_rank_in_category"])
            for metric in ("adjusted_held_out_rank",
                           "adjusted_held_out_rank_in_category"):
                med = (r.dropna(subset=[metric])
                         .groupby("true_kinase_node_id", observed=True)[metric]
                         .median())
                a = med[[k for k in med.index if k in ever]].to_numpy()
                b = med[[k for k in med.index if k not in ever]].to_numpy()
                if len(a) >= 5 and len(b) >= 5:
                    u, pv = stats.mannwhitneyu(a, b, alternative="two-sided")
                    rb = 1.0 - 2.0 * u / (len(a) * len(b))
                    rows.append({
                        "source": "LOO {}".format(run.replace("_multi_kinase_spectral", "")),
                        "metric": metric,
                        "n_ever_negative": len(a), "median_ever": float(np.median(a)),
                        "n_never_negative": len(b), "median_never": float(np.median(b)),
                        "mannwhitney_p": float(pv), "rank_biserial": float(rb),
                    })

        comp = pd.DataFrame(rows)
        if comp.empty:
            add("  no comparable rank source found")
        else:
            add("  Lower rank = better. Positive rank-biserial means the")
            add("  ever-negative group has the LARGER (worse) values.")
            add("")
            d = comp.copy()
            for c in ("median_ever", "median_never", "rank_biserial"):
                d[c] = d[c].map("{:.3f}".format)
            d["mannwhitney_p"] = d["mannwhitney_p"].map(
                lambda v: "{:.2e}".format(v) if v < 1e-3 else "{:.4f}".format(v))
            add(d.to_string(index=False))
            comp.to_csv(OUT / "negative_sampling_audit.csv", index=False)
            add("")

        # distribution of never-negative kinases across the alphabet
        never_idx = [order.index(k) for k in never if k in order]
        if never_idx:
            add("  never-negative kinases occupy candidate-list positions")
            add("  {} to {}, i.e. essentially everything after the head slice.".format(
                min(never_idx), max(never_idx)))
        add("")

    add("=" * 78)
    add("WHAT THE CODE DOES, STATED PLAINLY")
    add("=" * 78)
    add("  For every training site, the negative examples are the first 50")
    add("  alphabetically-ordered candidate kinases that are not known to act on")
    add("  that site. The same small set of early-alphabet kinases is therefore")
    add("  used as a negative example at nearly every site, and the large")
    add("  majority of candidates are never presented to the classifier as a")
    add("  negative at all.")
    add("")
    add("  This is a defect in the training design, not in the evaluation: the")
    add("  ranking still scores all 415 candidates. Its consequence is that the")
    add("  decision boundary is fitted against an alphabetically determined")
    add("  subset of the candidate space.")

    text = "\n".join(L)
    (OUT / "negative_sampling_audit.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
