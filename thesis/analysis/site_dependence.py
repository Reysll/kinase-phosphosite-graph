"""
Item 1 -- site-dependence diagnostic (read-only).

QUESTION
--------
Does the model actually rank kinases differently at different sites, or does it
apply a near-fixed popularity ordering everywhere? If kinase identity alone
explains almost all the variance in rank_in_site, then the "prediction" is
mostly a global prior and the site contributes little.

METRICS, per (network, embedding), over the Phase 3 predict-all matrices
    eta_squared   one-way ANOVA effect size for rank_in_site ~ kinase,
                  i.e. SS_between(kinase) / SS_total. 1.0 = rank is fully
                  determined by which kinase it is; 0.0 = kinase identity
                  explains nothing.
    R2            identical quantity expressed as regression R^2 (same number).
    median per-kinase SD of rank across sites
    median per-SITE SD of rank across kinases (reference scale)
    pct of kinases with SD <= 1 rank across all sites

Inputs:
    outputs_prediction/inference_all/{generic,control,liver_fc}_{spectral,node2vec}/ranks.csv.gz

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.site_dependence
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

INFER = Path("outputs_prediction/inference_all")
OUT = Path("thesis/analysis/out")

NETWORKS = {"Generic": "generic", "Control": "control", "Liver FC": "liver_fc"}
EMBEDDINGS = ["spectral", "node2vec"]


def eta_squared(values: np.ndarray, groups: np.ndarray) -> float:
    """One-way ANOVA eta^2 = SS_between / SS_total, groups given as integer codes."""
    grand = values.mean()
    ss_total = float(((values - grand) ** 2).sum())
    if ss_total == 0:
        return float("nan")
    counts = np.bincount(groups)
    sums = np.bincount(groups, weights=values)
    nz = counts > 0
    means = sums[nz] / counts[nz]
    ss_between = float((counts[nz] * (means - grand) ** 2).sum())
    return ss_between / ss_total


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L: list[str] = ["ITEM 1 -- SITE-DEPENDENCE DIAGNOSTIC", "=" * 78, ""]
    add = L.append
    add("eta^2 of rank_in_site ~ kinase identity.")
    add("  1.0 = a kinase gets the same rank at every site (pure global prior)")
    add("  0.0 = kinase identity carries no information about rank")
    add("")

    rows = []
    for emb in EMBEDDINGS:
        for label, d in NETWORKS.items():
            p = INFER / f"{d}_{emb}" / "ranks.csv.gz"
            if not p.exists():
                add(f"[missing] {p}")
                continue
            df = pd.read_csv(p, usecols=["site_node_id", "kinase_node_id", "rank_in_site"])
            vals = df["rank_in_site"].to_numpy(dtype=float)
            kcodes = pd.factorize(df["kinase_node_id"], sort=False)[0]

            e2 = eta_squared(vals, kcodes)
            per_kin = df.groupby("kinase_node_id")["rank_in_site"].std(ddof=1)
            per_site = df.groupby("site_node_id")["rank_in_site"].std(ddof=1)

            rows.append({
                "embedding": emb,
                "network": label,
                "n_sites": int(df["site_node_id"].nunique()),
                "n_kinases": int(df["kinase_node_id"].nunique()),
                "n_pairs": int(len(df)),
                "eta_squared": round(e2, 4),
                "R2": round(e2, 4),
                "median_SD_rank_per_kinase": round(float(per_kin.median()), 2),
                "mean_SD_rank_per_kinase": round(float(per_kin.mean()), 2),
                "pct_kinases_SD_le_1": round(100.0 * float((per_kin <= 1.0).mean()), 1),
                "pct_kinases_SD_le_5": round(100.0 * float((per_kin <= 5.0).mean()), 1),
                "median_SD_rank_per_site": round(float(per_site.median()), 2),
            })

    if not rows:
        print("\n".join(L))
        return

    t = pd.DataFrame(rows)
    add("RESULTS")
    add("-" * 78)
    add(t.to_string(index=False))
    t.to_csv(OUT / "site_dependence.csv", index=False)
    add("")
    add("Reference: if every kinase held a fixed rank at all sites, eta^2 = 1.000")
    add("and median per-kinase SD = 0.00. The per-SITE SD column is the spread of")
    add("ranks within one site, which is ~constant by construction (ranks 1..K).")

    text = "\n".join(L)
    (OUT / "site_dependence.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
