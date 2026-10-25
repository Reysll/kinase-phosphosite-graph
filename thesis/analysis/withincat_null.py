"""
Within-category comparative kinase ranking vs the degree-preserving null
(read-only; consumes existing outputs only).

QUESTION
--------
Dr. Ayati prefers comparative kinase ranking with kinases ranked inside their
own substrate-count category (poor / average / rich) over PSP-based prediction
performance. Earlier work showed the degree-preserving null reproduces the
Control -> Liver FC rank shift when kinases are ranked GLOBALLY. This script
asks the same question for WITHIN-CATEGORY ranking, for every kinase rather
than the six-kinase focus set of null_casestudy.py:

  How many kinases have a Control -> Liver FC within-category rank shift that
  differs significantly from the Control -> Null shift?

METHOD
------
1. Load Phase 3 predict-all ranks for Generic, Control and Liver FC (spectral)
   plus the degree-preserving null (spectral, seed 1).
2. Assign each kinase to poor / average / rich using the SAME rule as
   src_prediction/data_prep/kinase_categories.py (POOR_MAX=5, RICH_MIN=20),
   counting unique PSP substrate sites per kinase from the full generic edge
   table.
3. Re-rank within category at each site. Within a site, rank_in_site is a
   strictly monotone function of the model score, so ordering the members of
   one category by rank_in_site ascending reproduces exactly what re-scoring
   within the category would give, with no float-tie ambiguity.
4. Per kinase, take the median within-category rank across the sites common to
   Control / Liver FC / Null, and form the two shifts:
       liverfc_shift = median(Liver FC) - median(Control)
       null_shift    = median(Null)     - median(Control)
   Negative = the kinase ranks better (closer to 1) under that network.
5. Per kinase, ks_2samp on the within-category rank distributions, for
   Control vs Liver FC, Control vs Null, and Liver FC vs Null. The third is the
   contrast that answers the question: it asks whether the real network's
   within-category ranking is distinguishable from the randomized one.
   BH-FDR is applied across all tested kinases within each contrast, and counts
   are broken down by category. A per-category BH variant is also reported,
   since the category sets are the units Dr. Ayati cares about.

Inputs (all read-only):
    outputs/edges.csv.gz
    outputs_prediction/inference_all/{generic,control,liver_fc}_spectral/ranks.csv.gz
    thesis/analysis/out/null_inference/nullinfer_spectral_seed1_ranks.csv.gz

Writes only into thesis/analysis/out/:
    withincat_null.txt
    withincat_null_per_kinase.csv
    withincat_null_ks.csv
    withincat_null_summary.csv

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.withincat_null
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("thesis/analysis/out")
INFER = Path("outputs_prediction/inference_all")
NULL_RANKS = OUT / "null_inference" / "nullinfer_spectral_seed1_ranks.csv.gz"
GENERIC_EDGES = Path("outputs/edges.csv.gz")

POOR_MAX = 5
RICH_MIN = 20
CATEGORIES = ("poor", "average", "rich")
MIN_SITES = 30
COLS = ["site_node_id", "kinase_node_id", "rank_in_site"]


def psp_substrate_counts() -> dict:
    """Same rule as src_prediction/data_prep/kinase_categories.py, read-only."""
    edges = pd.read_csv(GENERIC_EDGES, usecols=["source", "target", "relation"])
    psp = edges.loc[edges["relation"] == "phosphorylates"]
    return psp.groupby("source")["target"].nunique().to_dict()


def assign(kinase_node_ids, counts: dict) -> pd.Series:
    out = {}
    for kid in kinase_node_ids:
        n = counts.get(kid, 0)
        if n <= POOR_MAX:
            out[kid] = "poor"
        elif n <= RICH_MIN:
            out[kid] = "average"
        else:
            out[kid] = "rich"
    return pd.Series(out, name="category")


def load_ranks(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=COLS)
    df["rank_in_site"] = df["rank_in_site"].astype(np.int32)
    return df


def add_within_category_rank(df: pd.DataFrame, cat_map: pd.Series) -> pd.DataFrame:
    """rank_in_cat = position within (site, category) ordered by rank_in_site."""
    df = df.copy()
    df["category"] = df["kinase_node_id"].map(cat_map)
    df = df.dropna(subset=["category"])
    df = df.sort_values(["site_node_id", "category", "rank_in_site"], kind="stable")
    grp = df.groupby(["site_node_id", "category"], observed=True)
    df["rank_in_cat"] = grp.cumcount().astype(np.int32) + 1
    return df


def fmt(v: float) -> str:
    return "{:.2e}".format(v) if v < 1e-3 else "{:.4f}".format(v)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L = []
    add = L.append
    add("WITHIN-CATEGORY COMPARATIVE RANKING vs DEGREE-PRESERVING NULL")
    add("=" * 78)
    add("")
    add("Read-only. Consumes existing Phase 3 and null-inference outputs.")
    add("Within-category rank: 1 = best inside that kinase's own category.")
    add("Shift convention: negative = ranks BETTER than under Control.")
    add("")

    sources = {
        "Generic": INFER / "generic_spectral" / "ranks.csv.gz",
        "Control": INFER / "control_spectral" / "ranks.csv.gz",
        "Liver FC": INFER / "liver_fc_spectral" / "ranks.csv.gz",
        "Null": NULL_RANKS,
    }
    missing = ["{}: {}".format(k, v) for k, v in sources.items() if not v.exists()]
    if missing:
        add("MISSING INPUTS:")
        for m in missing:
            add("  " + m)
        (OUT / "withincat_null.txt").write_text("\n".join(L), encoding="utf-8")
        print("\n".join(L))
        return

    counts = psp_substrate_counts()
    raw = {k: load_ranks(v) for k, v in sources.items()}

    all_kin = sorted(set(raw["Control"]["kinase_node_id"]))
    cat_map = assign(all_kin, counts)

    add("KINASE CATEGORY ASSIGNMENT (poor <=5, average 6-20, rich >20 PSP sites)")
    add("-" * 78)
    sizes = cat_map.value_counts().reindex(list(CATEGORIES)).fillna(0).astype(int)
    for c in CATEGORIES:
        add("  {:8s}: {:3d} kinases scored in the inference pass".format(c, sizes[c]))
    add("  {:8s}: {:3d}".format("total", len(cat_map)))
    add("")

    data = {k: add_within_category_rank(d, cat_map) for k, d in raw.items()}
    for k, d in data.items():
        add("{:9s}: {:,} sites, {:,} pairs".format(
            k, d["site_node_id"].nunique(), len(d)))
    add("")

    common_ctx = set.intersection(
        *[set(data[k]["site_node_id"]) for k in ("Control", "Liver FC", "Null")]
    )
    common_all = set.intersection(*[set(d["site_node_id"]) for d in data.values()])
    add("Sites common to Control / Liver FC / Null : {:,}".format(len(common_ctx)))
    add("Sites common to all four                  : {:,}".format(len(common_all)))
    add("  Contrasts below use the Control/Liver FC/Null intersection.")
    add("  Generic covers a smaller, PSP-annotated site population and is")
    add("  reported for reference only, never inside a contrast.")
    add("")

    ctx = {k: d[d["site_node_id"].isin(common_ctx)] for k, d in data.items()}
    med = {}
    for net in ("Control", "Liver FC", "Null"):
        med[net] = ctx[net].groupby("kinase_node_id", observed=True)["rank_in_cat"].median()
    gen = data["Generic"]
    gen = gen[gen["site_node_id"].isin(common_all)]
    med["Generic"] = gen.groupby("kinase_node_id", observed=True)["rank_in_cat"].median()

    per = pd.DataFrame(med)
    per["category"] = cat_map
    per["kinase"] = [str(k).split(":")[-1] for k in per.index]
    per["psp_sites"] = [counts.get(k, 0) for k in per.index]
    per["liverfc_shift"] = per["Liver FC"] - per["Control"]
    per["null_shift"] = per["Null"] - per["Control"]
    per["shift_gap"] = per["liverfc_shift"] - per["null_shift"]
    per = per.reset_index().rename(columns={"index": "kinase_node_id"})

    add("PER-KINASE WITHIN-CATEGORY MEDIAN SHIFTS, BY CATEGORY")
    add("-" * 78)
    agg = per.groupby("category", observed=True).agg(
        n_kinases=("kinase", "size"),
        median_liverfc_shift=("liverfc_shift", "median"),
        median_null_shift=("null_shift", "median"),
        median_shift_gap=("shift_gap", "median"),
        n_lf_better=("liverfc_shift", lambda s: int((s < 0).sum())),
        n_lf_worse=("liverfc_shift", lambda s: int((s > 0).sum())),
        n_lf_same=("liverfc_shift", lambda s: int((s == 0).sum())),
        n_null_better=("null_shift", lambda s: int((s < 0).sum())),
        n_null_worse=("null_shift", lambda s: int((s > 0).sum())),
        n_null_same=("null_shift", lambda s: int((s == 0).sum())),
    ).reindex(list(CATEGORIES))
    add(agg.to_string())
    add("")

    add("KS TESTS ON WITHIN-CATEGORY RANK DISTRIBUTIONS")
    add("-" * 78)
    add("One test per kinase with >= {} sites in both arms.".format(MIN_SITES))
    add("BH-FDR applied across ALL tested kinases within each contrast.")
    add("")

    piv = {}
    for net in ("Control", "Liver FC", "Null"):
        piv[net] = ctx[net].pivot_table(
            index="kinase_node_id", columns="site_node_id",
            values="rank_in_cat", observed=True)

    contrasts = [("Control", "Liver FC"), ("Control", "Null"), ("Liver FC", "Null")]
    ks_frames = []
    for a, b in contrasts:
        recs = []
        A, B = piv[a], piv[b]
        for kid in A.index.intersection(B.index):
            xa = A.loc[kid].dropna().to_numpy()
            xb = B.loc[kid].dropna().to_numpy()
            if len(xa) < MIN_SITES or len(xb) < MIN_SITES:
                continue
            ks, p = stats.ks_2samp(xa, xb)
            recs.append({
                "contrast": "{} vs {}".format(a, b),
                "kinase_node_id": kid,
                "kinase": str(kid).split(":")[-1],
                "category": cat_map.get(kid),
                "n_a": len(xa), "n_b": len(xb),
                "median_a": float(np.median(xa)),
                "median_b": float(np.median(xb)),
                "median_shift": float(np.median(xb) - np.median(xa)),
                "ks_statistic": float(ks),
                "p_value": float(p),
            })
        t = pd.DataFrame(recs)
        if t.empty:
            continue
        t["adj_p"] = stats.false_discovery_control(t["p_value"].to_numpy(), method="bh")
        t["significant"] = t["adj_p"] < 0.05
        t["adj_p_within_category"] = np.nan
        for c in t["category"].dropna().unique():
            m = t["category"] == c
            t.loc[m, "adj_p_within_category"] = stats.false_discovery_control(
                t.loc[m, "p_value"].to_numpy(), method="bh")
        t["significant_within_category"] = t["adj_p_within_category"] < 0.05
        ks_frames.append(t)

        add("--- {} vs {} ---".format(a, b))
        tab = t.groupby("category", observed=True).agg(
            n_tested=("kinase", "size"),
            n_significant=("significant", "sum"),
            n_sig_percat_bh=("significant_within_category", "sum"),
            median_ks=("ks_statistic", "median"),
            median_shift=("median_shift", "median"),
        ).reindex(list(CATEGORIES))
        tab["pct_significant"] = (100.0 * tab["n_significant"] / tab["n_tested"]).round(1)
        add(tab.to_string())
        add("  TOTAL tested {}, significant (global BH) {} ({:.1f}%)".format(
            len(t), int(t["significant"].sum()), 100.0 * t["significant"].mean()))
        add("")

    if ks_frames:
        full = pd.concat(ks_frames, ignore_index=True)
        full.to_csv(OUT / "withincat_null_ks.csv", index=False)

        add("HEADLINE: KINASES WHOSE LIVER FC RANKING DIFFERS FROM THE NULL")
        add("-" * 78)
        lf = full[full["contrast"] == "Liver FC vs Null"]
        if not lf.empty:
            tab = lf.groupby("category", observed=True).agg(
                n_tested=("kinase", "size"),
                n_differ_from_null=("significant", "sum"),
                n_differ_percat_bh=("significant_within_category", "sum"),
            ).reindex(list(CATEGORIES))
            tab["pct"] = (100.0 * tab["n_differ_from_null"] / tab["n_tested"]).round(1)
            add(tab.to_string())
            add("")
            add("  Across all categories: {} of {} kinases ({:.1f}%) have a".format(
                int(lf["significant"].sum()), len(lf),
                100.0 * lf["significant"].mean()))
            add("  within-category rank distribution under the real Liver FC network")
            add("  that is distinguishable from the degree-preserving null at")
            add("  BH-FDR < 0.05.")
            add("")
            top = lf.sort_values("adj_p").head(25)[
                ["kinase", "category", "n_a", "median_a", "median_b",
                 "median_shift", "ks_statistic", "adj_p"]]
            add("  25 kinases most distinguishable from the null:")
            disp = top.copy()
            for c in ("median_a", "median_b", "median_shift", "ks_statistic"):
                disp[c] = disp[c].map("{:.2f}".format)
            disp["adj_p"] = disp["adj_p"].map(fmt)
            add(disp.to_string(index=False))
            add("")

            add("EFFECT SIZES BEHIND THE KS SIGNIFICANCE")
            add("-" * 78)
            add("Each kinase is tested on 3,228 sites, so a KS test detects")
            add("distributional differences far smaller than one rank position.")
            add("The magnitude of the median shift is reported separately here.")
            add("")
            sig = lf[lf["significant"]]
            add("Distribution of |median_shift| among the {} significant kinases:".format(
                len(sig)))
            vc = sig["median_shift"].abs().value_counts().sort_index()
            for v, n in vc.items():
                add("  |median_shift| = {:>4}  : {:3d} kinases".format(v, int(n)))
            add("")
            for thr in (1, 2, 3, 5, 10):
                n = int((sig["median_shift"].abs() >= thr).sum())
                add("  |median_shift| >= {:2d} ranks : {:3d} of {}".format(
                    thr, n, len(sig)))
            add("")
            add("  ks_statistic among significant: median {:.3f}, max {:.3f}".format(
                sig["ks_statistic"].median(), sig["ks_statistic"].max()))
            add("")
            add("AGREEMENT BETWEEN THE REAL AND NULL SHIFTS (all {} kinases)".format(
                len(per)))
            add("-" * 78)
            for col in ("liverfc_shift", "null_shift", "shift_gap"):
                a = per[col].abs()
                add("  {:14s} median {:.1f}  mean {:.2f}  max {:.1f}  nonzero {:3d}".format(
                    col, a.median(), a.mean(), a.max(), int((a > 0).sum())))
            add("")
            add("  Pearson r (liverfc_shift vs null_shift)  = {:.4f}".format(
                per["liverfc_shift"].corr(per["null_shift"])))
            add("  Spearman  (liverfc_shift vs null_shift)  = {:.4f}".format(
                per["liverfc_shift"].corr(per["null_shift"], method="spearman")))
            add("")

        summ = full.groupby(["contrast", "category"], observed=True).agg(
            n_tested=("kinase", "size"),
            n_significant=("significant", "sum"),
            n_sig_percat_bh=("significant_within_category", "sum"),
            median_ks=("ks_statistic", "median"),
            median_shift=("median_shift", "median"),
        ).reset_index()
        summ.to_csv(OUT / "withincat_null_summary.csv", index=False)

    per.to_csv(OUT / "withincat_null_per_kinase.csv", index=False)

    add("FILES WRITTEN")
    add("-" * 78)
    for f in ("withincat_null.txt", "withincat_null_per_kinase.csv",
              "withincat_null_ks.csv", "withincat_null_summary.csv"):
        add("  thesis/analysis/out/" + f)

    text = "\n".join(L)
    (OUT / "withincat_null.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
