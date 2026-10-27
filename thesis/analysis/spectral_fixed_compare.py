"""
Compare the null-space-fixed spectral runs against the originals (read-only).

Reports, in the order requested, results before interpretation:

  1. degree vs mean rank (Spearman) and eta-squared, fixed vs original
  2. median per-kinase rank SD, fixed vs original
  3. paired Wilcoxon with rank-biserial: fixed vs original, and Generic vs
     Liver FC under the fix
  4. popularity baseline vs the fixed models, overall and by rich/average/poor

The only difference between the two sides is `drop_null_space`. Frozen trials,
candidate set, negative sampler, dimensions and random_state are identical, so
every LOO contrast is paired on trial_index.

Statistics are computed with the same helpers as the earlier rounds
(`stat_tests.paired_test`, `randneg_compare.eta_sq_and_sd`) so the numbers are
comparable to what is already in thesis/analysis/out/.

Inputs (read-only):
    thesis/analysis/out/spectral_fixed/fixed_loo_{generic,liver_fc}_spectral_results.csv.gz
    thesis/analysis/out/spectral_fixed/fixed_infer_{control,liver_fc}_spectral_ranks.csv.gz
    outputs_prediction/{generic,liver}_multi_kinase_spectral/results.csv.gz
    outputs_prediction/inference_all/{control,liver_fc}_spectral/ranks.csv.gz

Writes only into thesis/analysis/out/spectral_fixed/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.spectral_fixed_compare
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from thesis.analysis.degree_rank import _fit
from thesis.analysis.randneg_compare import ORIG_LOO, eta_sq_and_sd, ptest, show
from thesis.analysis.spectral_nullspace import _embedding_graph
from thesis.analysis.stat_tests import KEY, RANK_COL, holm

PRED = Path("outputs_prediction")
FX = Path("thesis/analysis/out/spectral_fixed")
CAT_COL = "held_out_kinase_category"
INCAT = "adjusted_held_out_rank_in_category"
CATS = ("rich", "average", "poor")

ORIG_INFER = {
    "control": PRED / "inference_all" / "control_spectral" / "ranks.csv.gz",
    "liver_fc": PRED / "inference_all" / "liver_fc_spectral" / "ranks.csv.gz",
}
LOO_COLS = [KEY, "true_kinase_node_id", RANK_COL, INCAT, CAT_COL]


def degree_and_eta(path: Path, deg: dict, sub_counts: dict) -> dict:
    """eta^2, per-kinase rank SD, and degree vs mean rank, from one ranks file."""
    out = eta_sq_and_sd(path)
    d = pd.read_csv(path, usecols=["kinase_node_id", "rank_in_site"])
    mean_rank = d.groupby("kinase_node_id")["rank_in_site"].mean()
    y = mean_rank.to_numpy(float)
    degrees = np.array([deg.get(k, np.nan) for k in mean_rank.index], dtype=float)
    subs = np.array([sub_counts.get(k, 0) for k in mean_rank.index], dtype=float)
    fd = _fit(degrees, y)
    fs = _fit(subs, y)
    out.update({
        "degree_rho": fd["spearman_rho"], "degree_r2": fd["pearson_r2"],
        "substrates_rho": fs["spearman_rho"], "substrates_r2": fs["pearson_r2"],
    })
    return out


def main() -> None:
    from src_prediction.core.config import CANDIDATE_KINASES_OUT, GENERIC_EDGES
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.data_prep.kinase_categories import (
        assign_kinase_categories, compute_psp_substrate_counts)

    FX.mkdir(parents=True, exist_ok=True)
    L: list[str] = []
    add = L.append
    add("SPECTRAL NULL-SPACE FIX vs ORIGINAL SPECTRAL RUNS")
    add("=" * 78)
    add("")
    add("Read-only. `drop_null_space=True` is the only difference: same frozen")
    add("trials, same 415 scored candidates, same original negative sampler,")
    add("spectral d=32 random_state=42. Every LOO contrast is paired on")
    add("trial_index. Results first, interpretation last.")
    add("")

    cands = read_csv_auto(CANDIDATE_KINASES_OUT)
    counts = compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES))
    cats = assign_kinase_categories(cands["node_id"].astype(str).tolist(), counts)

    # ---------------------------------------------- 1 + 2. predict-all side
    add("1 + 2. DEGREE vs MEAN RANK, ETA-SQUARED AND PER-KINASE RANK SD")
    add("-" * 78)
    add("From the predict-all output. degree_rho is Spearman between a kinase's")
    add("embedding-graph degree and its mean rank_in_site over all scored sites.")
    add("eta^2 = 1 means a kinase holds the same rank at every site.")
    add("")
    rows1 = []
    for net in ("control", "liver_fc"):
        fixed = FX / "fixed_infer_{}_spectral_ranks.csv.gz".format(net)
        if not fixed.exists():
            add("  [missing] {}".format(fixed))
            continue
        deg = dict(_embedding_graph(net).degree())
        for run, p in (("original", ORIG_INFER[net]), ("fixed", fixed)):
            if not p.exists():
                add("  [missing] {}".format(p))
                continue
            r = degree_and_eta(p, deg, counts)
            r.update({"network": net, "run": run})
            rows1.append(r)
            print("[cmp] {:<9} {:<8}: degree_rho={:+.4f} eta2={:.4f} sd={:.2f}".format(
                net, run, r["degree_rho"], r["eta_squared"],
                r["median_SD_rank_per_kinase"]))
    if rows1:
        d1 = pd.DataFrame(rows1)
        pd.DataFrame(rows1).to_csv(FX / "fixed_degree_eta.csv", index=False)
        disp = d1[["network", "run", "n_sites", "n_kinases", "degree_rho", "degree_r2",
                   "substrates_rho", "eta_squared", "median_SD_rank_per_kinase",
                   "pct_kinases_SD_le_1"]].copy()
        for c in ("degree_rho", "substrates_rho"):
            disp[c] = disp[c].map("{:+.4f}".format)
        for c in ("degree_r2", "eta_squared"):
            disp[c] = disp[c].map("{:.4f}".format)
        for c in ("median_SD_rank_per_kinase", "pct_kinases_SD_le_1"):
            disp[c] = disp[c].map("{:.2f}".format)
        add(disp.to_string(index=False))
    add("")
    add("  Reference, original spectral: degree_rho -0.9995 (control) and")
    add("  -0.9973 (liver FC), eta^2 = 0.9997, median per-kinase SD 1.05 and")
    add("  1.18. source: thesis/analysis/out/relational_only/degree_rank.txt,")
    add("  thesis/analysis/out/site_dependence.txt")
    add("")

    # ---------------------------------------------- load LOO frames
    new_loo, orig_loo = {}, {}
    for net in ("generic", "liver_fc"):
        p = FX / "fixed_loo_{}_spectral_results.csv.gz".format(net)
        if p.exists():
            new_loo[net] = pd.read_csv(p, usecols=LOO_COLS)
            orig_loo[net] = pd.read_csv(ORIG_LOO[net], usecols=LOO_COLS)
    if not new_loo:
        add("NO FIXED LOO RUNS FOUND under {}".format(FX))
        txt = "\n".join(L)
        (FX / "spectral_fixed_compare.txt").write_text(txt, encoding="utf-8")
        print(txt)
        return

    bad: set[int] = set()
    for df in list(new_loo.values()) + list(orig_loo.values()):
        bad |= set(df.loc[df[RANK_COL].isna(), KEY].astype(int))
    keep = set(orig_loo[list(orig_loo)[0]][KEY].astype(int)) - bad
    add("3. PAIRED TESTS ON THE LOO RANKS")
    add("-" * 78)
    add("Trials rankable in every run compared here: {:,}".format(len(keep)))
    add("")

    # ---------------------------------------------- 3a. fixed vs original
    add("3a. ORIGINAL vs FIXED, per network")
    add("Convention: diff = rank_A - rank_B; positive favours B = fixed.")
    add("")
    rows3 = []
    for net, nd in sorted(new_loo.items()):
        od = orig_loo[net]
        m = (od[od[KEY].isin(keep)][[KEY, RANK_COL, INCAT, CAT_COL]]
             .rename(columns={RANK_COL: "a_glob", INCAT: "a_cat"})
             .merge(nd[nd[KEY].isin(keep)][[KEY, RANK_COL, INCAT]]
                    .rename(columns={RANK_COL: "b_glob", INCAT: "b_cat"}), on=KEY))
        for ca, cb, metric in (("a_glob", "b_glob", "global"),
                               ("a_cat", "b_cat", "within-category")):
            sub = m.dropna(subset=[ca, cb])
            r = ptest(sub[ca].to_numpy(float), sub[cb].to_numpy(float),
                      "original", "fixed")
            r.update({"network": net, "metric": metric, "stratum": "ALL", "n": len(sub)})
            rows3.append(r)
            for c in CATS:
                s2 = sub[sub[CAT_COL] == c]
                if len(s2) < 10:
                    continue
                r2 = ptest(s2[ca].to_numpy(float), s2[cb].to_numpy(float),
                           "original", "fixed")
                r2.update({"network": net, "metric": metric, "stratum": c, "n": len(s2)})
                rows3.append(r2)
    ps = holm([r["p_value"] for r in rows3])
    for r, p in zip(rows3, ps):
        r["p_holm"] = p
    cols = [c for c in ["network", "metric", "stratum", "n", "median_A", "median_B",
                        "median_paired_diff", "rank_biserial", "B_wins", "A_wins",
                        "n_ties", "p_value", "p_holm", "favours"] if c in rows3[0]]
    add(show(rows3, cols))
    pd.DataFrame(rows3).to_csv(FX / "fixed_vs_original.csv", index=False)
    add("")

    # ---------------------------------------------- 3b. context effect
    add("3b. CONTEXT EFFECT UNDER THE FIX (Generic vs Liver FC)")
    add("")
    rows4 = []
    if "generic" in new_loo and "liver_fc" in new_loo:
        for label, src in (("fixed", new_loo), ("original (reference)", orig_loo)):
            g, lv = src["generic"], src["liver_fc"]
            m = (g[g[KEY].isin(keep)][[KEY, RANK_COL, CAT_COL]].rename(columns={RANK_COL: "a"})
                 .merge(lv[lv[KEY].isin(keep)][[KEY, RANK_COL]].rename(columns={RANK_COL: "b"}),
                        on=KEY)).dropna()
            r = ptest(m["a"].to_numpy(float), m["b"].to_numpy(float), "Generic", "Liver FC")
            r.update({"run": label, "stratum": "ALL", "n": len(m)})
            rows4.append(r)
            for c in CATS:
                s2 = m[m[CAT_COL] == c]
                if len(s2) < 10:
                    continue
                r2 = ptest(s2["a"].to_numpy(float), s2["b"].to_numpy(float),
                           "Generic", "Liver FC")
                r2.update({"run": label, "stratum": c, "n": len(s2)})
                rows4.append(r2)
        c4 = [c for c in ["run", "stratum", "n", "median_A", "median_B",
                          "median_paired_diff", "rank_biserial", "B_wins", "A_wins",
                          "n_ties", "p_value", "favours"] if c in rows4[0]]
        add(show(rows4, c4))
        pd.DataFrame(rows4).to_csv(FX / "fixed_context_effect.csv", index=False)
    else:
        add("  needs both generic and liver_fc fixed LOO runs")
    add("")
    add("  Reference: the original Generic vs Liver FC spectral contrast was")
    add("  median 79.0 -> 74.0, rank-biserial +0.380, p_holm = 1.40e-162.")
    add("  source: thesis/analysis/out/stat_tests.txt")
    add("")

    # ---------------------------------------------- 4. popularity baseline
    add("4. POPULARITY BASELINE vs THE FIXED MODELS")
    add("-" * 78)
    add("Baseline ranks kinases by substrate count over the 6,909 training")
    add("edges (model parity), no embedding. Convention: positive diff favours")
    add("B = the fixed model.")
    add("")
    try:
        from thesis.analysis.popularity_baseline import build_baseline, scored_candidate_set

        base = build_baseline(keep, cats, counts_from="frozen",
                              candidates=scored_candidate_set())
        rows5 = []
        for net, nd in sorted(new_loo.items()):
            m = (base[base[KEY].isin(keep)][[KEY, RANK_COL, "rank_in_category", "category"]]
                 .rename(columns={RANK_COL: "pop_glob", "rank_in_category": "pop_cat"})
                 .merge(nd[nd[KEY].isin(keep)][[KEY, RANK_COL, INCAT]]
                        .rename(columns={RANK_COL: "mod_glob", INCAT: "mod_cat"}), on=KEY))
            for pa, pb, metric in (("pop_glob", "mod_glob", "global"),
                                   ("pop_cat", "mod_cat", "within-category")):
                sub = m.dropna(subset=[pa, pb])
                r = ptest(sub[pa].to_numpy(float), sub[pb].to_numpy(float),
                          "popularity", "fixed model")
                r.update({"network": net, "metric": metric, "stratum": "ALL",
                          "n": len(sub)})
                rows5.append(r)
                for c in CATS:
                    s2 = sub[sub["category"] == c]
                    if len(s2) < 10:
                        continue
                    r2 = ptest(s2[pa].to_numpy(float), s2[pb].to_numpy(float),
                               "popularity", "fixed model")
                    r2.update({"network": net, "metric": metric, "stratum": c,
                               "n": len(s2)})
                    rows5.append(r2)
        c5 = [c for c in ["network", "metric", "stratum", "n", "median_A", "median_B",
                          "median_paired_diff", "rank_biserial", "B_wins", "A_wins",
                          "n_ties", "p_value", "favours"] if c in rows5[0]]
        add(show(rows5, c5))
        pd.DataFrame(rows5).to_csv(FX / "fixed_vs_popularity.csv", index=False)
        add("")
        add("  Reference, original spectral: popularity median 25.0 against")
        add("  Generic 79.0 (rank-biserial -0.758) and Liver FC 74.0 (-0.686).")
        add("  source: thesis/analysis/out/popularity_baseline.txt")
    except Exception as e:                                        # noqa: BLE001
        add("  baseline comparison FAILED: {}".format(e))
    add("")

    add("FILES WRITTEN")
    add("-" * 78)
    for f in ("spectral_fixed_compare.txt", "fixed_degree_eta.csv",
              "fixed_vs_original.csv", "fixed_context_effect.csv",
              "fixed_vs_popularity.csv"):
        add("  thesis/analysis/out/spectral_fixed/" + f)

    txt = "\n".join(L)
    (FX / "spectral_fixed_compare.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
