"""
Verification suite requested after Checkpoint 2.

READ-ONLY. Reads outputs_prediction/ and imports project modules for category
definitions only. Writes solely into thesis/analysis/out/.

Steps
-----
V1. Kinase-level tests using MEDIAN adjusted rank per kinase (vs MEAN previously).
V2. Stratified paired tests within `rich` and `poor` kinase categories
    (definitions per category_ranks.txt / kinase_categories.py), Holm-corrected
    within each stratum, at trial level and kinase level.
V3. Informative-restricted (G15) analysis at kinase level.
V4. Sanity outputs: per-kinase node2vec-vs-spectral scatter, improved/worsened
    counts per comparison, histogram of trials per kinase.
V5. Integrity check: identical trial sets, identical kinase sets, identical
    missing-rank handling across all eight runs.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.verify_tests

Seed: none required (deterministic).
"""

from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from thesis.analysis.stat_tests import (
    EMBEDDINGS,
    KEY,
    NETWORKS,
    RANK_COL,
    build_family,
    holm,
    informative_kinases,
    load_all,
    paired_test,
    union_unrankable,
)

OUT = Path("thesis/analysis/out")
PRED = Path("outputs_prediction")


# --------------------------------------------------------------------- helpers
def categories() -> dict[str, str]:
    """kinase_node_id -> poor/average/rich, using the project's own definitions."""
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.core.config import CANDIDATE_KINASES_OUT, GENERIC_EDGES
    from src_prediction.data_prep.kinase_categories import (
        assign_kinase_categories,
        compute_psp_substrate_counts,
    )

    generic_edges = read_csv_auto(GENERIC_EDGES)
    counts = compute_psp_substrate_counts(generic_edges)
    cands = read_csv_auto(CANDIDATE_KINASES_OUT)["node_id"].astype(str).tolist()
    return assign_kinase_categories(cands, counts)


def fmt(df: pd.DataFrame) -> str:
    d = df.copy()
    for c in ("median_A", "median_B", "median_paired_diff", "wilcoxon_W"):
        if c in d:
            d[c] = d[c].map(lambda v: f"{v:,.1f}")
    if "rank_biserial" in d:
        d["rank_biserial"] = d["rank_biserial"].map(lambda v: f"{v:+.3f}")
    for c in ("p_value", "p_holm"):
        if c in d:
            d[c] = d[c].map(lambda v: f"{v:.3e}" if v < 1e-3 else f"{v:.4f}")
    return d.to_string(index=False)


def run_family(
    frames, keep, *, level: str, agg: str = "mean",
    kinase_filter: set[str] | None = None,
) -> pd.DataFrame:
    rows = []
    for block, na, ea, nb, eb in build_family():
        da, db = frames[(na, ea)], frames[(nb, eb)]
        da = da[da[KEY].isin(keep)]
        db = db[db[KEY].isin(keep)]
        if kinase_filter is not None:
            da = da[da["true_kinase_node_id"].isin(kinase_filter)]
            db = db[db["true_kinase_node_id"].isin(kinase_filter)]
        if da.empty or db.empty:
            continue
        if level == "kinase":
            ga = da.groupby("true_kinase_node_id")[RANK_COL].agg(agg)
            gb = db.groupby("true_kinase_node_id")[RANK_COL].agg(agg)
            common = ga.index.intersection(gb.index)
            if len(common) < 3:
                continue
            a, b = ga.loc[common].to_numpy(), gb.loc[common].to_numpy()
        else:
            m = da.merge(db, on=KEY, suffixes=("_a", "_b"))
            a, b = m[f"{RANK_COL}_a"].to_numpy(), m[f"{RANK_COL}_b"].to_numpy()
        res = paired_test(a, b, f"{na} ({ea})", f"{nb} ({eb})")
        res["block"] = block
        # improved / worsened counts (lower rank = improved for B)
        diff = a - b
        res["n_B_better"] = int((diff > 0).sum())
        res["n_A_better"] = int((diff < 0).sum())
        rows.append(res)
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["p_holm"] = holm(out["p_value"].tolist())
    out["significant_holm_0.05"] = out["p_holm"] < 0.05
    return out[[
        "block", "comparison", "n", "n_ties", "median_A", "median_B",
        "median_paired_diff", "rank_biserial", "n_B_better", "n_A_better",
        "p_value", "p_holm", "significant_holm_0.05", "favours",
    ]]


# --------------------------------------------------------------------- main
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames, prov = load_all()
    bad = union_unrankable(frames)
    all_trials = set(frames[("Generic", "node2vec")][KEY].astype(int))
    keep = all_trials - bad
    cats = categories()

    L: list[str] = []
    add = L.append
    add("VERIFICATION SUITE -- post-Checkpoint-2 requests")
    add("=" * 78)
    add("")

    # ---------------------------------------------------------------- V5 first
    add("V5. INTEGRITY CHECK ACROSS THE EIGHT RUNS")
    add("-" * 78)
    ref_key = ("Generic", "node2vec")
    ref = frames[ref_key].sort_values(KEY).reset_index(drop=True)
    trial_ok, kinase_ok, nan_ok = True, True, True
    nan_sets = {}
    for k, df in frames.items():
        d = df.sort_values(KEY).reset_index(drop=True)
        trial_ok &= set(d[KEY]) == set(ref[KEY])
        kinase_ok &= d["true_kinase_node_id"].equals(ref["true_kinase_node_id"])
        nan_sets[k] = set(d.loc[d[RANK_COL].isna(), KEY].astype(int))
    all_nan_same = len({frozenset(v) for v in nan_sets.values()}) == 1
    add(f"  Identical trial_index set in all 8 runs      : {trial_ok}")
    add(f"  Identical held-out kinase per trial_index    : {kinase_ok}")
    add(f"  Identical NaN(adjusted rank) trial set in all: {all_nan_same}")
    add(f"  NaN count per run                            : "
        f"{sorted({len(v) for v in nan_sets.values()})}")
    add(f"  Union of unrankable trials                   : {len(bad)}")
    add(f"  Trials retained for every test               : {len(keep):,}")
    kset = {k: set(df[df[KEY].isin(keep)]['true_kinase_node_id']) for k, df in frames.items()}
    add(f"  Identical kinase set after filtering         : "
        f"{len({frozenset(v) for v in kset.values()}) == 1}")
    add(f"  Distinct held-out kinases                    : {len(kset[ref_key])}")
    add("")

    # ---------------------------------------------------------------- V1
    add("V1. KINASE-LEVEL USING MEDIAN (previously MEAN)")
    add("-" * 78)
    v1_med = run_family(frames, keep, level="kinase", agg="median")
    v1_mean = run_family(frames, keep, level="kinase", agg="mean")
    add("MEDIAN aggregation:")
    add(fmt(v1_med))
    add("")
    add("MEAN aggregation (for direct comparison):")
    add(fmt(v1_mean[["comparison", "median_paired_diff", "rank_biserial",
                     "p_holm", "significant_holm_0.05", "favours"]]))
    add("")
    chg = []
    for c in v1_med["comparison"]:
        a = v1_med.loc[v1_med["comparison"] == c].iloc[0]
        b = v1_mean.loc[v1_mean["comparison"] == c].iloc[0]
        if bool(a["significant_holm_0.05"]) != bool(b["significant_holm_0.05"]) or \
           a["favours"] != b["favours"]:
            chg.append(f"    {c}: median->({a['favours']}, sig={a['significant_holm_0.05']})"
                       f" vs mean->({b['favours']}, sig={b['significant_holm_0.05']})")
    add("  Conclusions changed vs MEAN? " + ("YES" if chg else "NO -- identical direction and significance"))
    for c in chg:
        add(c)
    add("")

    # ---------------------------------------------------------------- V2
    add("V2. STRATIFIED BY KINASE CATEGORY (Holm within each stratum)")
    add("-" * 78)
    ncat = pd.Series(cats).value_counts().to_dict()
    add(f"  Kinase counts by category (all 420 candidates): {ncat}")
    for stratum in ("rich", "poor"):
        ks = {k for k, c in cats.items() if c == stratum}
        n_tr = int(frames[ref_key][frames[ref_key][KEY].isin(keep)]
                   ["true_kinase_node_id"].isin(ks).sum())
        add("")
        add(f"  --- {stratum.upper()} (kinases={len(ks)}, trials={n_tr:,}) --- trial level")
        t = run_family(frames, keep, level="trial", kinase_filter=ks)
        add(fmt(t) if not t.empty else "    (insufficient data)")
        if not t.empty:
            t.to_csv(OUT / f"verify_stratified_{stratum}_trial.csv", index=False)
        add(f"  --- {stratum.upper()} --- kinase level (median agg)")
        tk = run_family(frames, keep, level="kinase", agg="median", kinase_filter=ks)
        add(fmt(tk) if not tk.empty else "    (insufficient data)")
        if not tk.empty:
            tk.to_csv(OUT / f"verify_stratified_{stratum}_kinase.csv", index=False)
    add("")

    # ---------------------------------------------------------------- V3
    add("V3. INFORMATIVE-RESTRICTED (G15) AT KINASE LEVEL (median agg)")
    add("-" * 78)
    inf_n2v, inf_spec = informative_kinases("node2vec"), informative_kinases("spectral")
    inf_both = inf_n2v & inf_spec
    add(f"  informative(node2vec)={len(inf_n2v)}  informative(spectral)={len(inf_spec)}"
        f"  intersection={len(inf_both)}")
    rows = []
    for block, na, ea, nb, eb in build_family():
        filt = inf_both if block.startswith("B") else informative_kinases(ea)
        sub = run_family(frames, keep, level="kinase", agg="median", kinase_filter=filt)
        sel = sub[sub["comparison"] == f"{na} ({ea}) vs {nb} ({eb})"]
        if not sel.empty:
            rows.append(sel.iloc[0])
    v3 = pd.DataFrame(rows).reset_index(drop=True)
    v3["p_holm"] = holm(v3["p_value"].tolist())
    v3["significant_holm_0.05"] = v3["p_holm"] < 0.05
    add(fmt(v3))
    v3.to_csv(OUT / "verify_informative_kinase_level.csv", index=False)
    add("")

    # ---------------------------------------------------------------- V4
    add("V4. SANITY OUTPUTS")
    add("-" * 78)
    add("  Improved / worsened kinase counts (kinase level, median agg):")
    add(fmt(v1_med[["comparison", "n", "n_B_better", "n_A_better", "n_ties", "favours"]]))
    add("")
    tr_per_kin = (frames[ref_key][frames[ref_key][KEY].isin(keep)]
                  .groupby("true_kinase_node_id").size().sort_values(ascending=False))
    add(f"  Trials per kinase: n={len(tr_per_kin)}  min={tr_per_kin.min()}  "
        f"median={int(tr_per_kin.median())}  mean={tr_per_kin.mean():.1f}  max={tr_per_kin.max()}")
    add(f"  Top 5 by trial count: "
        f"{', '.join(f'{k.split(':')[-1]}={v}' for k, v in tr_per_kin.head(5).items())}")
    share = tr_per_kin.sort_values(ascending=False).cumsum() / tr_per_kin.sum()
    add(f"  Kinases covering 50% of all trials: {int((share <= 0.5).sum()) + 1}")
    add(f"  Kinases with a single trial        : {int((tr_per_kin == 1).sum())}")
    tr_per_kin.rename("n_trials").to_csv(OUT / "verify_trials_per_kinase.csv")

    # scatter + histogram
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.4), dpi=150)
    net = "Liver FC"
    ga = (frames[(net, "node2vec")][frames[(net, "node2vec")][KEY].isin(keep)]
          .groupby("true_kinase_node_id")[RANK_COL].median())
    gb = (frames[(net, "spectral")][frames[(net, "spectral")][KEY].isin(keep)]
          .groupby("true_kinase_node_id")[RANK_COL].median())
    common = ga.index.intersection(gb.index)
    cc = [cats.get(k, "poor") for k in common]
    colors = {"poor": "#d62728", "average": "#ff7f0e", "rich": "#2b6cb0"}
    for cat in ("rich", "average", "poor"):
        m = [i for i, c in enumerate(cc) if c == cat]
        axes[0].scatter(ga.loc[common].to_numpy()[m], gb.loc[common].to_numpy()[m],
                        s=18, alpha=0.7, c=colors[cat], label=f"{cat} (n={len(m)})")
    axes[0].plot([0, 420], [0, 420], "k--", lw=1, alpha=0.6)
    axes[0].set_xlabel("node2vec median adjusted rank"); axes[0].set_ylabel("spectral median adjusted rank")
    axes[0].set_title(f"Per-kinase rank: node2vec vs spectral ({net})\nbelow diagonal = spectral better")
    axes[0].legend(fontsize=9); axes[0].set_xlim(0, 420); axes[0].set_ylim(0, 420)

    axes[1].hist(tr_per_kin.to_numpy(), bins=60, color="#2b6cb0", edgecolor="white")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("trials per kinase"); axes[1].set_ylabel("number of kinases (log)")
    axes[1].set_title(f"Trial-count imbalance across kinases\n(n={len(tr_per_kin)}, max={tr_per_kin.max()})")
    plt.tight_layout()
    plt.savefig(OUT / "verify_sanity_plots.png", bbox_inches="tight")
    plt.close()
    add(f"  Wrote scatter + histogram: {OUT / 'verify_sanity_plots.png'}")
    add("")

    v1_med.to_csv(OUT / "verify_kinase_median.csv", index=False)
    text = "\n".join(L)
    (OUT / "verify_tests.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
