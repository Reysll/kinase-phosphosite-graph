"""
Item 3 -- case-study robustness against the degree-preserving null (read-only).

QUESTION
--------
The necroptosis cluster (RIPK1, RIPK3, MLKL, IRAK4, TRPM7) and TTK were
identified from Phase 3 predict-all rank shifts between Control and Liver FC
under spectral embedding. If the SAME shifts appear when the Liver FC
correlation edges are replaced by degree-matched random edges, the shift is a
property of the added degree structure rather than of liver co-variation.

Compares per-kinase rank_in_site distributions across the 3,228 FC-valid sites:
    Generic, Control, real Liver FC, Null(Liver FC randomized)

Tests, per kinase (scipy.stats.ks_2samp, BH-FDR across the kinase set):
    Control  vs Liver FC   -- the original published contrast
    Control  vs Null       -- does the null reproduce it?
    Liver FC vs Null       -- do real and null differ from each other?

Inputs:
    outputs_prediction/inference_all/{generic,control,liver_fc}_spectral/ranks.csv.gz
    thesis/analysis/out/null_inference/nullinfer_spectral_seed1_ranks.csv.gz

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.null_casestudy
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("thesis/analysis/out")
INFER = Path("outputs_prediction/inference_all")
NULL_RANKS = Path("thesis/analysis/out/null_inference/nullinfer_spectral_seed1_ranks.csv.gz")

FOCUS = ["TTK", "RIPK1", "RIPK3", "MLKL", "IRAK4", "TRPM7"]
COLS = ["site_node_id", "kinase_node_id", "rank_in_site"]


def load(path: Path, label: str) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=COLS)
    df["kinase"] = df["kinase_node_id"].astype(str).str.split(":").str[-1]
    df["network"] = label
    return df


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    sources = {
        "Generic": INFER / "generic_spectral" / "ranks.csv.gz",
        "Control": INFER / "control_spectral" / "ranks.csv.gz",
        "Liver FC": INFER / "liver_fc_spectral" / "ranks.csv.gz",
        "Null": NULL_RANKS,
    }
    L: list[str] = ["ITEM 3 -- CASE-STUDY ROBUSTNESS vs DEGREE-PRESERVING NULL",
                    "=" * 78, ""]
    add = L.append

    missing = [f"{k}: {v}" for k, v in sources.items() if not v.exists()]
    if missing:
        add("MISSING INPUTS:")
        for m in missing:
            add("  " + m)
        print("\n".join(L))
        return

    data = {k: load(v, k) for k, v in sources.items()}
    for k, d in data.items():
        add(f"{k:9s}: {d['site_node_id'].nunique():,} sites, {len(d):,} pairs")
    add("")

    # restrict to sites present in all four, so distributions are comparable
    site_sets = [set(d["site_node_id"]) for d in data.values()]
    common_all = set.intersection(*site_sets)
    common_ctx = set.intersection(*[set(data[k]["site_node_id"])
                                    for k in ("Control", "Liver FC", "Null")])
    add(f"Sites common to all four networks        : {len(common_all):,}")
    add(f"Sites common to Control/Liver FC/Null    : {len(common_ctx):,}")
    add("  (Generic covers only its own smaller site population; contrasts below")
    add("   use the Control/Liver FC/Null intersection, Generic shown for reference)")
    add("")

    # ---------------------------------------------------------- median table
    add("MEDIAN rank_in_site PER KINASE (lower = predicted more relevant)")
    add("-" * 78)
    rows = []
    for kin in FOCUS:
        r = {"kinase": kin}
        for net, d in data.items():
            sub = d[(d["kinase"] == kin)]
            sites = common_all if net == "Generic" else common_ctx
            sub = sub[sub["site_node_id"].isin(sites)]
            r[net] = float(sub["rank_in_site"].median()) if len(sub) else np.nan
            if net == "Null":
                r["n_sites"] = int(len(sub))
        r["LiverFC_minus_Control"] = r["Liver FC"] - r["Control"]
        r["Null_minus_Control"] = r["Null"] - r["Control"]
        rows.append(r)
    med = pd.DataFrame(rows)[
        ["kinase", "n_sites", "Generic", "Control", "Liver FC", "Null",
         "LiverFC_minus_Control", "Null_minus_Control"]
    ]
    add(med.to_string(index=False))
    med.to_csv(OUT / "null_casestudy_medians.csv", index=False)
    add("")

    # ---------------------------------------------------------- KS tests
    add("KS-TESTS PER KINASE (BH-FDR across the 6 focus kinases, per contrast)")
    add("-" * 78)
    contrasts = [("Control", "Liver FC"), ("Control", "Null"), ("Liver FC", "Null")]
    out_rows = []
    for a, b in contrasts:
        recs = []
        for kin in FOCUS:
            xa = data[a][(data[a]["kinase"] == kin) &
                         (data[a]["site_node_id"].isin(common_ctx))]["rank_in_site"]
            xb = data[b][(data[b]["kinase"] == kin) &
                         (data[b]["site_node_id"].isin(common_ctx))]["rank_in_site"]
            if len(xa) < 30 or len(xb) < 30:
                continue
            ks, p = stats.ks_2samp(xa.to_numpy(), xb.to_numpy())
            recs.append({"contrast": f"{a} vs {b}", "kinase": kin,
                         "n_a": len(xa), "n_b": len(xb),
                         "median_a": float(xa.median()), "median_b": float(xb.median()),
                         "median_shift": float(xb.median() - xa.median()),
                         "ks_statistic": float(ks), "p_value": float(p),
                         "std_a": float(xa.std()), "std_b": float(xb.std())})
        if not recs:
            continue
        t = pd.DataFrame(recs)
        t["adj_p"] = stats.false_discovery_control(t["p_value"].to_numpy(), method="bh")
        t["significant"] = t["adj_p"] < 0.05
        out_rows.append(t)
        disp = t.copy()
        for c in ("median_a", "median_b", "median_shift", "ks_statistic", "std_a", "std_b"):
            disp[c] = disp[c].map("{:.2f}".format)
        for c in ("p_value", "adj_p"):
            disp[c] = disp[c].map(lambda v: f"{v:.2e}" if v < 1e-3 else f"{v:.4f}")
        add(disp.to_string(index=False))
        add("")

    if out_rows:
        full = pd.concat(out_rows, ignore_index=True)
        full.to_csv(OUT / "null_casestudy_ks.csv", index=False)

        add("DIRECT QUESTION: does the null reproduce the Liver FC shift?")
        add("-" * 78)
        piv = full.pivot_table(index="kinase", columns="contrast",
                               values="median_shift")
        if {"Control vs Liver FC", "Control vs Null"} <= set(piv.columns):
            piv["null_reproduces_pct"] = (
                100.0 * piv["Control vs Null"] / piv["Control vs Liver FC"]
            ).round(1)
            add(piv.round(2).to_string())
            add("")
            add("  null_reproduces_pct = (Null shift) / (Liver FC shift) x 100.")
            add("  ~100 means the randomized network produces the same shift;")
            add("  ~0 means the shift requires the real correlation structure.")

    text = "\n".join(L)
    (OUT / "null_casestudy.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
