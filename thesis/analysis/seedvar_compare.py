"""
Seed-variance comparison for node2vec (read-only).

PURPOSE
-------
node2vec is stochastic. Every headline node2vec number in the project comes from
a single seed (42). This script measures how far the metric moves when ONLY the
seed changes, and puts that seed-to-seed effect size directly beside the two
substantive node2vec effect sizes:

    Generic vs Liver FC   (the context effect, RQ1)
    Null    vs Liver FC   (the degree-preserving null contrast, Tier 2)

If seed-to-seed movement is comparable to either, that contrast cannot be
interpreted without replication. If it is much smaller, the contrast stands.

INPUTS (must already exist; this script does not run any experiment)
    thesis/analysis/out/seed_variance/seedvar_{net}_node2vec_seed{S}_results.csv.gz
        produced by: thesis.analysis.seed_variance  (step b of run_overnight.ps1)
    outputs_prediction/{generic,liver}_multi_kinase/results.csv.gz   (seed 42)
    thesis/analysis/out/null/null_node2vec_seed*_results.csv.gz      (optional)

METRICS, per (network, seed-vs-seed42) pair
    * Spearman rho of per-trial adjusted_held_out_rank
    * paired Wilcoxon + matched-pairs rank-biserial + win/loss/tie counts
    * top-1 prediction agreement rate across trials
    * overlap of the top-10 kinase leaderboard (kinases ordered by mean
      adjusted rank), which is the LOO-side analogue of the Phase 3
      top-10-per-site overlap used in tier1_analyses.py

NOTE ON SETTINGS PARITY
    seed_variance.py pins every Node2VecStrategy argument to the values in
    run_multi_kinase_{generic,liver}.py (notably batch_words=10000) so that the
    seed is the only thing that differs. Runs made before that fix are
    quarantined and must not be fed to this script.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.seedvar_compare
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from thesis.analysis.stat_tests import KEY, RANK_COL, paired_test
from thesis.analysis.verify_tests import fmt

SEED_DIR = Path("thesis/analysis/out/seed_variance")
NULL_DIR = Path("thesis/analysis/out/null")
PRED = Path("outputs_prediction")
OUT = Path("thesis/analysis/out")

BASELINE = {  # the original seed-42 node2vec runs
    "generic": PRED / "generic_multi_kinase" / "results.csv.gz",
    "liver_fc": PRED / "liver_multi_kinase" / "results.csv.gz",
}
TOP_N = 10


def load(path: Path) -> pd.DataFrame:
    cols = [KEY, "true_kinase_node_id", RANK_COL]
    df = pd.read_csv(path)
    if "top1_predicted_kinase" in df.columns:
        cols.append("top1_predicted_kinase")
    return df[cols]


def leaderboard(df: pd.DataFrame, n: int = TOP_N) -> list[str]:
    """Top-n kinases by mean adjusted rank (lower = better)."""
    return (df.groupby("true_kinase_node_id")[RANK_COL]
              .mean().sort_values().head(n).index.tolist())


def compare(a: pd.DataFrame, b: pd.DataFrame, la: str, lb: str) -> dict:
    m = a.merge(b, on=KEY, suffixes=("_a", "_b")).dropna(
        subset=[f"{RANK_COL}_a", f"{RANK_COL}_b"]
    )
    x = m[f"{RANK_COL}_a"].to_numpy()
    y = m[f"{RANK_COL}_b"].to_numpy()
    res = paired_test(x, y, la, lb)
    d = x - y
    res["B_wins"] = int((d > 0).sum())
    res["A_wins"] = int((d < 0).sum())
    res["ties"] = int((d == 0).sum())
    res["pct_ties"] = round(100.0 * (d == 0).mean(), 2)
    rho, _ = stats.spearmanr(x, y)
    res["spearman_rho"] = float(rho)

    if {"top1_predicted_kinase_a", "top1_predicted_kinase_b"} <= set(m.columns):
        agree = (m["top1_predicted_kinase_a"] == m["top1_predicted_kinase_b"]).mean()
        res["top1_agreement_pct"] = round(100.0 * float(agree), 2)
    else:
        res["top1_agreement_pct"] = float("nan")

    la_top, lb_top = leaderboard(a), leaderboard(b)
    res["top10_leaderboard_overlap"] = len(set(la_top) & set(lb_top))
    return res


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L: list[str] = ["NODE2VEC SEED-VARIANCE vs SUBSTANTIVE EFFECTS", "=" * 78, ""]
    add = L.append
    add("Convention: diff = rank_A - rank_B; positive favours B.")
    add("Seed-to-seed rows are NOISE estimates, not findings.")
    add("")

    rows: list[dict] = []

    # ---- seed-to-seed, within network ------------------------------------
    found_any = False
    for net, base_path in BASELINE.items():
        if not base_path.exists():
            add(f"[missing baseline] {base_path}")
            continue
        base = load(base_path)
        for p in sorted(SEED_DIR.glob(f"seedvar_{net}_node2vec_seed*_results.csv.gz")):
            s = re.search(r"seed(\d+)_results", p.name)
            if not s:
                continue
            found_any = True
            r = compare(base, load(p), f"{net} seed42", f"{net} seed{s.group(1)}")
            r["family"] = "seed-to-seed (noise)"
            rows.append(r)

    if not found_any:
        add("No seed-variance runs found. Run step b of run_overnight.ps1 first:")
        add("  powershell -ExecutionPolicy Bypass -File thesis\\analysis\\run_overnight.ps1 -Steps b")
        add("")

    # ---- substantive contrasts -------------------------------------------
    if all(p.exists() for p in BASELINE.values()):
        gen, liv = load(BASELINE["generic"]), load(BASELINE["liver_fc"])
        r = compare(gen, liv, "Generic seed42", "Liver FC seed42")
        r["family"] = "context effect (RQ1)"
        rows.append(r)

        for p in sorted(NULL_DIR.glob("null_node2vec_seed*_results.csv.gz")):
            s = re.search(r"seed(\d+)_results", p.name)
            r = compare(load(p), liv, f"Null seed{s.group(1)}", "Liver FC seed42")
            r["family"] = "null contrast (Tier 2)"
            rows.append(r)

    if not rows:
        print("\n".join(L))
        return

    t = pd.DataFrame(rows)
    cols = ["family", "comparison", "n", "median_A", "median_B",
            "median_paired_diff", "rank_biserial", "spearman_rho",
            "top1_agreement_pct", "top10_leaderboard_overlap",
            "B_wins", "A_wins", "ties", "pct_ties", "p_value"]
    t = t[[c for c in cols if c in t.columns]]
    t = t.sort_values("family")

    add("ALL CONTRASTS SIDE BY SIDE")
    add("-" * 78)
    add(fmt(t))
    add("")

    noise = t[t["family"].str.startswith("seed-to-seed")]
    if not noise.empty:
        add("INTERPRETATION GUIDE (numbers only; no claim asserted here)")
        add("-" * 78)
        add(f"  max |rank_biserial| from seed alone      : "
            f"{noise['rank_biserial'].abs().max():.3f}")
        add(f"  max |median paired diff| from seed alone : "
            f"{noise['median_paired_diff'].abs().max():.1f}")
        add("  Compare each substantive contrast above against these ceilings.")

    text = "\n".join(L)
    (OUT / "seedvar_compare.txt").write_text(text, encoding="utf-8")
    t.to_csv(OUT / "seedvar_compare.csv", index=False)
    print(text)


if __name__ == "__main__":
    main()
