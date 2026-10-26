"""
Compare the random-negatives runs against the originals (read-only).

Answers, in the order requested:
  1. Original vs random negatives, paired Wilcoxon with rank-biserial and
     win/loss/tie, per network, overall and by rich/average/poor.
  2. Generic vs Liver FC under random negatives: does the context effect change?
  3. eta-squared and per-kinase rank SD from the new predict-all output, against
     the original 0.9997.
  4. The popularity baseline against the random-negative models, overall and
     stratified.
  5. How many distinct kinases appear as negatives under the new sampler.

Inputs (read-only):
    thesis/analysis/out/random_negatives/randneg_*_results.csv.gz
    thesis/analysis/out/random_negatives/randneg_infer_*_ranks.csv.gz
    thesis/analysis/out/random_negatives/randneg_*_negative_usage.csv
    outputs_prediction/{generic,liver}_multi_kinase_spectral/results.csv.gz
    outputs_prediction/inference_all/{...}_spectral/ranks.csv.gz

Writes only into thesis/analysis/out/random_negatives/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.randneg_compare
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from thesis.analysis.stat_tests import KEY, RANK_COL, holm, paired_test

PRED = Path("outputs_prediction")
RN = Path("thesis/analysis/out/random_negatives")
CAT_COL = "held_out_kinase_category"
INCAT = "adjusted_held_out_rank_in_category"
CATS = ("rich", "average", "poor")

# original spectral LOO runs, keyed by the network name used in the new runs
ORIG_LOO = {
    "generic": PRED / "generic_multi_kinase_spectral" / "results.csv.gz",
    "liver_fc": PRED / "liver_multi_kinase_spectral" / "results.csv.gz",
    "control": PRED / "control_multi_kinase_spectral" / "results.csv.gz",
}
ORIG_INFER = {
    "generic": PRED / "inference_all" / "generic_spectral" / "ranks.csv.gz",
    "control": PRED / "inference_all" / "control_spectral" / "ranks.csv.gz",
    "liver_fc": PRED / "inference_all" / "liver_fc_spectral" / "ranks.csv.gz",
}


def eta_sq_and_sd(path: Path) -> dict:
    d = pd.read_csv(path, usecols=["site_node_id", "kinase_node_id", "rank_in_site"])
    y = d["rank_in_site"].to_numpy(float)
    grand = y.mean()
    g = d.groupby("kinase_node_id", observed=True)["rank_in_site"]
    means = g.mean()
    counts = g.size()
    ss_between = float((counts * (means - grand) ** 2).sum())
    ss_total = float(((y - grand) ** 2).sum())
    sd = g.std()
    return {
        "n_sites": int(d["site_node_id"].nunique()),
        "n_kinases": int(d["kinase_node_id"].nunique()),
        "n_pairs": int(len(d)),
        "eta_squared": ss_between / ss_total if ss_total else np.nan,
        "median_SD_rank_per_kinase": float(sd.median()),
        "mean_SD_rank_per_kinase": float(sd.mean()),
        "pct_kinases_SD_le_1": float(100 * (sd <= 1).mean()),
        "pct_kinases_SD_le_5": float(100 * (sd <= 5).mean()),
    }


def ptest(a, b, la, lb) -> dict:
    """paired_test plus explicit win/loss counts (lower rank wins)."""
    r = paired_test(a, b, la, lb)
    diff = np.asarray(a, float) - np.asarray(b, float)
    r["B_wins"] = int((diff > 0).sum())
    r["A_wins"] = int((diff < 0).sum())
    return r


def fmt_p(v: float) -> str:
    return "{:.3e}".format(v) if v < 1e-3 else "{:.4f}".format(v)


def show(rows: list[dict], cols: list[str]) -> str:
    d = pd.DataFrame(rows)[cols].copy()
    for c in ("p_value", "p_holm"):
        if c in d.columns:
            d[c] = d[c].map(fmt_p)
    for c in ("rank_biserial", "median_paired_diff", "median_A", "median_B"):
        if c in d.columns:
            d[c] = d[c].map(lambda v: "{:+.3f}".format(v) if isinstance(v, float) else v)
    return d.to_string(index=False)


def main() -> None:
    L = []
    add = L.append
    add("RANDOM NEGATIVES vs ORIGINAL HEAD-SLICE NEGATIVES")
    add("=" * 78)
    add("")
    add("Read-only comparison. The runs hold the frozen trial set, the candidate")
    add("set, the embeddings and the classifier fixed; only negative selection")
    add("differs, so every contrast below is paired on trial_index.")
    add("")

    # ------------------------------------------------- 5. negative coverage
    add("5. NEGATIVE COVERAGE UNDER THE NEW SAMPLER")
    add("-" * 78)
    real = set(pd.read_csv(PRED / "candidate_kinases.csv.gz")["node_id"].astype(str))
    cov = []
    for f in sorted(RN.glob("*_negative_usage.csv")):
        u = pd.read_csv(f)
        u.columns = ["kinase", "n_uses"]
        keep_u = u[u["kinase"].isin(real)]
        cov.append({"run": f.name.replace("_negative_usage.csv", ""),
                    "distinct_kinases": len(keep_u),
                    "total_pairs": int(keep_u["n_uses"].sum()),
                    "min_uses": int(keep_u["n_uses"].min()),
                    "median_uses": int(keep_u["n_uses"].median()),
                    "max_uses": int(keep_u["n_uses"].max()),
                    "logged_rows_excluded": int(len(u) - len(keep_u))})
    if cov:
        add(pd.DataFrame(cov).to_string(index=False))
    add("")
    add("  logged_rows_excluded counts synthetic ids written by an earlier")
    add("  version of verify_patch(); they never entered training. Counts above")
    add("  are restricted to the 420 real candidates.")
    add("")
    add("  Original head-slice sampler: 55 distinct kinases of 420,")
    add("  361 of 415 scored kinases (87.0%) never used as a negative.")
    add("  source: thesis/analysis/out/negative_sampling_audit.txt")
    add("")

    # ------------------------------------------------- load LOO frames
    new_loo = {}
    for f in sorted(RN.glob("randneg_loo_*_results.csv.gz")):
        parts = f.name.split("_")
        net = "liver_fc" if "liver" in f.name else parts[2]
        seed = int(f.name.split("seed")[1].split("_")[0])
        new_loo[(net, seed)] = pd.read_csv(
            f, usecols=[KEY, "true_kinase_node_id", RANK_COL, INCAT, CAT_COL])

    if not new_loo:
        add("NO LOO RUNS FOUND under {}".format(RN))
        (RN / "randneg_compare.txt").write_text("\n".join(L), encoding="utf-8")
        print("\n".join(L))
        return

    orig = {}
    for net in {k[0] for k in new_loo}:
        orig[net] = pd.read_csv(
            ORIG_LOO[net], usecols=[KEY, "true_kinase_node_id", RANK_COL, INCAT, CAT_COL])

    # shared rankable trials across everything loaded
    bad = set()
    for df in list(new_loo.values()) + list(orig.values()):
        bad |= set(df.loc[df[RANK_COL].isna(), KEY].astype(int))
    keep = set(orig[list(orig)[0]][KEY].astype(int)) - bad
    add("Trials rankable in every run compared here: {:,}".format(len(keep)))
    add("")

    # ------------------------------- 1. original vs random, per network
    add("1. ORIGINAL vs RANDOM NEGATIVES (paired, per network)")
    add("-" * 78)
    add("Convention: diff = rank_A - rank_B; positive favours B = random negatives.")
    add("")
    rows = []
    for (net, seed), nd in sorted(new_loo.items()):
        od = orig[net]
        m = (od[od[KEY].isin(keep)][[KEY, RANK_COL, INCAT, CAT_COL, "true_kinase_node_id"]]
             .rename(columns={RANK_COL: "a_glob", INCAT: "a_cat"})
             .merge(nd[nd[KEY].isin(keep)][[KEY, RANK_COL, INCAT]]
                    .rename(columns={RANK_COL: "b_glob", INCAT: "b_cat"}), on=KEY))
        for col_a, col_b, metric in (("a_glob", "b_glob", "global"),
                                     ("a_cat", "b_cat", "within-category")):
            sub = m.dropna(subset=[col_a, col_b])
            r = ptest(sub[col_a].to_numpy(float), sub[col_b].to_numpy(float),
                            "original", "random-neg")
            r.update({"network": net, "seed": seed, "metric": metric,
                      "stratum": "ALL", "n": len(sub)})
            rows.append(r)
            for c in CATS:
                s2 = sub[sub[CAT_COL] == c]
                if len(s2) < 10:
                    continue
                r2 = ptest(s2[col_a].to_numpy(float), s2[col_b].to_numpy(float),
                                 "original", "random-neg")
                r2.update({"network": net, "seed": seed, "metric": metric,
                           "stratum": c, "n": len(s2)})
                rows.append(r2)
    ps = holm([r["p_value"] for r in rows])
    for r, p in zip(rows, ps):
        r["p_holm"] = p
    cols = ["network", "seed", "metric", "stratum", "n", "median_A", "median_B",
            "median_paired_diff", "rank_biserial", "B_wins", "A_wins", "n_ties",
            "p_value", "p_holm", "favours"]
    cols = [c for c in cols if c in rows[0]]
    add(show(rows, cols))
    pd.DataFrame(rows).to_csv(RN / "randneg_vs_original.csv", index=False)
    add("")

    # ------------------------------- 2. context effect under random negatives
    add("2. CONTEXT EFFECT UNDER RANDOM NEGATIVES (Generic vs Liver FC)")
    add("-" * 78)
    rows2 = []
    seeds = sorted({s for (_, s) in new_loo})
    for seed in seeds:
        if ("generic", seed) not in new_loo or ("liver_fc", seed) not in new_loo:
            continue
        g = new_loo[("generic", seed)]
        lv = new_loo[("liver_fc", seed)]
        m = (g[g[KEY].isin(keep)][[KEY, RANK_COL, CAT_COL]].rename(columns={RANK_COL: "a"})
             .merge(lv[lv[KEY].isin(keep)][[KEY, RANK_COL]].rename(columns={RANK_COL: "b"}),
                    on=KEY)).dropna()
        r = ptest(m["a"].to_numpy(float), m["b"].to_numpy(float),
                        "Generic(rand)", "Liver FC(rand)")
        r.update({"seed": seed, "stratum": "ALL", "n": len(m)})
        rows2.append(r)
        for c in CATS:
            s2 = m[m[CAT_COL] == c]
            if len(s2) < 10:
                continue
            r2 = ptest(s2["a"].to_numpy(float), s2["b"].to_numpy(float),
                             "Generic(rand)", "Liver FC(rand)")
            r2.update({"seed": seed, "stratum": c, "n": len(s2)})
            rows2.append(r2)
    # original comparison for reference
    if "generic" in orig and "liver_fc" in orig:
        m = (orig["generic"][orig["generic"][KEY].isin(keep)][[KEY, RANK_COL, CAT_COL]]
             .rename(columns={RANK_COL: "a"})
             .merge(orig["liver_fc"][orig["liver_fc"][KEY].isin(keep)][[KEY, RANK_COL]]
                    .rename(columns={RANK_COL: "b"}), on=KEY)).dropna()
        r = ptest(m["a"].to_numpy(float), m["b"].to_numpy(float),
                        "Generic(orig)", "Liver FC(orig)")
        r.update({"seed": 0, "stratum": "ALL (ORIGINAL, reference)", "n": len(m)})
        rows2.append(r)
    if rows2:
        c2 = [c for c in ["seed", "stratum", "n", "median_A", "median_B",
                          "median_paired_diff", "rank_biserial", "B_wins", "A_wins",
                          "n_ties", "p_value", "favours"] if c in rows2[0]]
        add(show(rows2, c2))
        pd.DataFrame(rows2).to_csv(RN / "randneg_context_effect.csv", index=False)
    add("")
    add("  The original Generic vs Liver FC spectral contrast was")
    add("  median 79.0 -> 74.0, rank-biserial +0.380, p_holm = 1.40e-162.")
    add("  source: thesis/analysis/out/stat_tests.txt")
    add("")

    # ------------------------------- 3. eta squared from new predict-all
    add("3. SITE DEPENDENCE FROM THE NEW PREDICT-ALL OUTPUT")
    add("-" * 78)
    rows3 = []
    for f in sorted(RN.glob("randneg_infer_*_ranks.csv.gz")):
        net = "liver_fc" if "liver" in f.name else f.name.split("_")[2]
        s = eta_sq_and_sd(f)
        s.update({"run": "random-neg", "network": net})
        rows3.append(s)
        if net in ORIG_INFER and ORIG_INFER[net].exists():
            s0 = eta_sq_and_sd(ORIG_INFER[net])
            s0.update({"run": "original", "network": net})
            rows3.append(s0)
    if rows3:
        d3 = pd.DataFrame(rows3)[
            ["network", "run", "n_sites", "n_kinases", "eta_squared",
             "median_SD_rank_per_kinase", "mean_SD_rank_per_kinase",
             "pct_kinases_SD_le_1", "pct_kinases_SD_le_5"]].sort_values(["network", "run"])
        for c in ("eta_squared",):
            d3[c] = d3[c].map("{:.4f}".format)
        for c in ("median_SD_rank_per_kinase", "mean_SD_rank_per_kinase",
                  "pct_kinases_SD_le_1", "pct_kinases_SD_le_5"):
            d3[c] = d3[c].map("{:.2f}".format)
        add(d3.to_string(index=False))
        pd.DataFrame(rows3).to_csv(RN / "randneg_site_dependence.csv", index=False)
    else:
        add("  no predict-all runs found yet")
    add("")
    add("  Original spectral Control and Liver FC: eta^2 = 0.9997,")
    add("  median per-kinase rank SD 1.05 and 1.18.")
    add("  source: thesis/analysis/out/site_dependence.txt")
    add("")

    # ------------------------------- 4. popularity baseline
    add("4. POPULARITY BASELINE vs RANDOM-NEGATIVE MODELS")
    add("-" * 78)
    try:
        from thesis.analysis.popularity_baseline import build_baseline, scored_candidate_set
        from src_prediction.core.io_utils import read_csv_auto
        from src_prediction.core.config import GENERIC_EDGES, CANDIDATE_KINASES_OUT
        from src_prediction.data_prep.kinase_categories import (
            assign_kinase_categories, compute_psp_substrate_counts)

        cands = read_csv_auto(CANDIDATE_KINASES_OUT)
        counts = compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES))
        cats = assign_kinase_categories(cands["node_id"].astype(str).tolist(), counts)
        base = build_baseline(keep, cats, counts_from="frozen",
                              candidates=scored_candidate_set())
        rows4 = []
        for (net, seed), nd in sorted(new_loo.items()):
            m = (base[base[KEY].isin(keep)][[KEY, RANK_COL, "rank_in_category", "category"]]
                 .rename(columns={RANK_COL: "pop_glob", "rank_in_category": "pop_cat"})
                 .merge(nd[nd[KEY].isin(keep)][[KEY, RANK_COL, INCAT]]
                        .rename(columns={RANK_COL: "mod_glob", INCAT: "mod_cat"}), on=KEY))
            for pa, pb, metric in (("pop_glob", "mod_glob", "global"),
                                   ("pop_cat", "mod_cat", "within-category")):
                sub = m.dropna(subset=[pa, pb])
                r = ptest(sub[pa].to_numpy(float), sub[pb].to_numpy(float),
                                "popularity", "random-neg model")
                r.update({"network": net, "seed": seed, "metric": metric,
                          "stratum": "ALL", "n": len(sub)})
                rows4.append(r)
                for c in CATS:
                    s2 = sub[sub["category"] == c]
                    if len(s2) < 10:
                        continue
                    r2 = ptest(s2[pa].to_numpy(float), s2[pb].to_numpy(float),
                                     "popularity", "random-neg model")
                    r2.update({"network": net, "seed": seed, "metric": metric,
                               "stratum": c, "n": len(s2)})
                    rows4.append(r2)
        if rows4:
            c4 = [c for c in ["network", "seed", "metric", "stratum", "n",
                              "median_A", "median_B", "median_paired_diff",
                              "rank_biserial", "B_wins", "A_wins", "n_ties",
                              "p_value", "favours"] if c in rows4[0]]
            add(show(rows4, c4))
            pd.DataFrame(rows4).to_csv(RN / "randneg_vs_popularity.csv", index=False)
        add("")
        add("  Original: popularity median 25.0 against Generic spectral 79.0")
        add("  (rank-biserial -0.758) and Liver FC spectral 74.0 (-0.686).")
        add("  source: thesis/analysis/out/popularity_baseline.txt")
    except Exception as e:                                        # noqa: BLE001
        add("  baseline comparison FAILED: {}".format(e))
    add("")

    add("FILES WRITTEN")
    add("-" * 78)
    for f in ("randneg_compare.txt", "randneg_vs_original.csv",
              "randneg_context_effect.csv", "randneg_site_dependence.csv",
              "randneg_vs_popularity.csv"):
        add("  thesis/analysis/out/random_negatives/" + f)

    text = "\n".join(L)
    (RN / "randneg_compare.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
