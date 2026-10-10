"""
Tier 2 analysis -- compare the degree-preserving null against real Liver FC and Generic.

READ-ONLY. Reads thesis/analysis/out/null/*_results.csv.gz (produced by
null_network.py) plus the project's existing LOO results. Writes only into
thesis/analysis/out/.

Comparisons, per embedding and per seed, on matched trials:
    Generic   vs Null        -- does adding degree-matched RANDOM edges help at all?
    Null      vs Liver FC    -- does the measured co-variation structure add
                                anything beyond edge count and degree?
    Generic   vs Liver FC    -- reference (reproduces the headline test)

Paired Wilcoxon, Holm-corrected within each (embedding, seed) family of 3,
with median paired difference and matched-pairs rank-biserial correlation.
Trials unrankable in ANY of the compared runs are dropped, and the same trial
set is used across the three tests of a family.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.null_compare

Seed: none required (deterministic given the inputs).
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pandas as pd

from thesis.analysis.stat_tests import KEY, RANK_COL, holm, paired_test
from thesis.analysis.verify_tests import categories, fmt

NULL_DIR = Path("thesis/analysis/out/null")
PRED = Path("outputs_prediction")
OUT = Path("thesis/analysis/out")

REAL = {
    ("spectral", "Generic"): PRED / "generic_multi_kinase_spectral" / "results.csv.gz",
    ("spectral", "Liver FC"): PRED / "liver_multi_kinase_spectral" / "results.csv.gz",
    ("node2vec", "Generic"): PRED / "generic_multi_kinase" / "results.csv.gz",
    ("node2vec", "Liver FC"): PRED / "liver_multi_kinase" / "results.csv.gz",
}


def load_rank(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)[[KEY, "true_kinase_node_id", RANK_COL]]


def discover_nulls() -> dict[str, list[tuple[int, Path]]]:
    found: dict[str, list[tuple[int, Path]]] = {}
    for p in sorted(NULL_DIR.glob("null_*_results.csv.gz")):
        m = re.match(r"null_(\w+?)_seed(\d+)_results\.csv\.gz", p.name)
        if not m or "probe" in p.name:
            continue
        emb, seed = m.group(1), int(m.group(2))
        found.setdefault(emb, []).append((seed, p))
    return found


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    nulls = discover_nulls()
    if not nulls:
        print("No completed null runs found in", NULL_DIR)
        return

    cats = categories()
    L: list[str] = ["TIER 2 -- DEGREE-PRESERVING NULL NETWORK", "=" * 78, ""]
    add = L.append
    add("Null = Liver FC correlation edges randomized by double-edge swap.")
    add("Degree sequence preserved exactly; ~99.9% of pairings destroyed.")
    add("Convention: positive median diff / positive r = B ranks BETTER.")
    add("")

    all_rows = []
    for emb, items in sorted(nulls.items()):
        gen = load_rank(REAL[(emb, "Generic")])
        liv = load_rank(REAL[(emb, "Liver FC")])
        for seed, p in sorted(items):
            nul = load_rank(p)
            if len(nul) < 6000:
                add(f"--- {emb} seed {seed}: SKIPPED (only {len(nul)} trials, run incomplete)")
                continue
            bad = set()
            for d in (gen, liv, nul):
                bad |= set(d.loc[d[RANK_COL].isna(), KEY].astype(int))
            keep = set(gen[KEY].astype(int)) - bad

            def arr(d: pd.DataFrame) -> np.ndarray:
                return (d[d[KEY].isin(keep)].sort_values(KEY)[RANK_COL].to_numpy())

            g, lv, nl = arr(gen), arr(liv), arr(nul)
            pairs = [
                (g, nl, "Generic", f"Null(s{seed})"),
                (nl, lv, f"Null(s{seed})", "Liver FC"),
                (g, lv, "Generic", "Liver FC"),
            ]
            rows = []
            for a, b, la, lb in pairs:
                r = paired_test(a, b, la, lb)
                d = a - b
                # win/loss/tie on the raw per-trial ranks (lower rank = win for B)
                r["B_wins"] = int((d > 0).sum())
                r["A_wins"] = int((d < 0).sum())
                r["ties"] = int((d == 0).sum())
                r["pct_ties"] = round(100.0 * (d == 0).mean(), 2)
                r["pct_B_wins"] = round(100.0 * (d > 0).mean(), 2)
                rows.append(r)
            t = pd.DataFrame(rows)
            t["p_holm"] = holm(t["p_value"].tolist())
            t["significant_holm_0.05"] = t["p_holm"] < 0.05
            t.insert(0, "seed", seed)
            t.insert(0, "embedding", emb)
            disp = t[["embedding", "seed", "comparison", "n", "median_A", "median_B",
                      "median_paired_diff", "rank_biserial", "B_wins", "A_wins",
                      "ties", "pct_ties", "p_holm", "significant_holm_0.05",
                      "favours"]]
            add(f"--- {emb} | seed {seed} | n = {len(g):,} matched trials ---")
            add(fmt(disp))
            add("")
            all_rows.append(t)

            # per-category medians for this seed
            kin = gen[gen[KEY].isin(keep)].sort_values(KEY)["true_kinase_node_id"].to_numpy()
            cat = np.array([cats.get(k, "poor") for k in kin])
            crows = []
            for c in ("rich", "average", "poor"):
                m = cat == c
                crows.append({"category": c, "n": int(m.sum()),
                              "Generic": float(np.median(g[m])),
                              "Null": float(np.median(nl[m])),
                              "Liver FC": float(np.median(lv[m]))})
            add("  median adjusted rank by category:")
            add("  " + pd.DataFrame(crows).to_string(index=False).replace("\n", "\n  "))
            add("")

    if all_rows:
        full = pd.concat(all_rows, ignore_index=True)
        full.to_csv(OUT / "null_compare.csv", index=False)
        add("SUMMARY -- 'Null vs Liver FC' across seeds")
        add("-" * 78)
        sub = full[full["comparison"].str.contains("vs Liver FC") &
                   full["comparison"].str.startswith("Null")]
        if not sub.empty:
            add(fmt(sub[["embedding", "seed", "n", "median_A", "median_B",
                         "median_paired_diff", "rank_biserial", "B_wins", "A_wins",
                         "ties", "pct_ties", "pct_B_wins", "p_holm",
                         "significant_holm_0.05", "favours"]]))
    # provenance
    meta = []
    for p in sorted(NULL_DIR.glob("null_*_meta.json")):
        if "probe" not in p.name:
            meta.append(json.loads(p.read_text(encoding="utf-8")))
    if meta:
        add("")
        add("NULL RUN METADATA")
        add("-" * 78)
        for m in meta:
            add(f"  {m['tag']}: median={m['adj_median_rank']}, mean={m['adj_mean_rank']:.1f}, "
                f"n={m['n_trials']}, unrankable={m['n_unrankable']}, "
                f"loo={m['loo_seconds']}s")

    text = "\n".join(L)
    (OUT / "null_compare.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
