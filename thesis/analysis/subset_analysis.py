"""
Item 7 -- evaluation-subset analysis (read-only, existing results only).

QUESTION
--------
The full 6,909-trial evaluation set is PhosphoSitePlus-derived and is not
restricted to anything measured in the liver proteomics. Liver context might
therefore be diluted: most trials concern sites and kinases the liver dataset
never observed. This restricts the same contrasts to trials where liver
measurement actually exists.

Subsets
    (a) SITE-measured : held-out site is among the 3,228 FC-valid liver sites
    (b) KINASE-measured: the held-out true kinase is detected in the liver
        proteomics (protein expression or phosphorylation sheet)
    (c) BOTH

Contrasts (each Holm-corrected within its own subset family)
    spectral : Generic  vs Liver FC
    spectral : Null(s1) vs Liver FC        [+ seeds 2,3 if present]
    node2vec : Generic  vs Liver FC on SEED-AVERAGED ranks {42,7,13}

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.subset_analysis
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

from thesis.analysis.stat_tests import KEY, RANK_COL, holm, paired_test
from thesis.analysis.verify_tests import fmt

PRED = Path("outputs_prediction")
NULL_DIR = Path("thesis/analysis/out/null")
SEED_DIR = Path("thesis/analysis/out/seed_variance")
OUT = Path("thesis/analysis/out")
SEEDS = [42, 7, 13]


def load(path: Path) -> pd.DataFrame:
    return pd.read_csv(path)[[KEY, "true_kinase_node_id", "site_node_id", RANK_COL]]


def liver_measured_sites() -> set[str]:
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values
    return set(load_site_fc_values(LIVER_EXCEL).index.astype(str))


def liver_measured_proteins() -> set[str]:
    """Gene symbols present in the liver graph as liver-observed protein nodes."""
    from src_prediction.core.config import LIVER_NODES
    from src_prediction.core.io_utils import read_csv_auto
    nodes = read_csv_auto(LIVER_NODES)
    prot = nodes.loc[nodes["node_type"] == "protein", "node_id"].astype(str)
    return set(prot)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L: list[str] = ["ITEM 7 -- EVALUATION-SUBSET ANALYSIS", "=" * 78, ""]
    add = L.append

    gen_s = load(PRED / "generic_multi_kinase_spectral" / "results.csv.gz")
    liv_s = load(PRED / "liver_multi_kinase_spectral" / "results.csv.gz")
    nulls = {}
    for p in sorted(NULL_DIR.glob("null_spectral_seed*_results.csv.gz")):
        s = re.search(r"seed(\d+)_results", p.name)
        if s and "probe" not in p.name:
            nulls[int(s.group(1))] = load(p)

    n2v: dict[str, dict[int, pd.DataFrame]] = {"generic": {}, "liver_fc": {}}
    base = {"generic": PRED / "generic_multi_kinase" / "results.csv.gz",
            "liver_fc": PRED / "liver_multi_kinase" / "results.csv.gz"}
    for net in n2v:
        for s in SEEDS:
            p = base[net] if s == 42 else SEED_DIR / f"seedvar_{net}_node2vec_seed{s}_results.csv.gz"
            if p.exists():
                n2v[net][s] = load(p)

    # drop trials unrankable anywhere
    allframes = [gen_s, liv_s, *nulls.values()] + [d for v in n2v.values() for d in v.values()]
    bad: set[int] = set()
    for d in allframes:
        bad |= set(d.loc[d[RANK_COL].isna(), KEY].astype(int))
    keep_all = set(gen_s[KEY].astype(int)) - bad

    ref = gen_s[gen_s[KEY].isin(keep_all)].sort_values(KEY).reset_index(drop=True)
    fc_sites = liver_measured_sites()
    liv_prot = liver_measured_proteins()

    site_ok = ref["site_node_id"].astype(str).isin(fc_sites).to_numpy()
    kin_ok = ref["true_kinase_node_id"].astype(str).isin(liv_prot).to_numpy()
    trial_ids = ref[KEY].to_numpy()

    subsets = {
        "FULL (all rankable)": set(trial_ids),
        "(a) site is FC-valid": set(trial_ids[site_ok]),
        "(b) kinase liver-detected": set(trial_ids[kin_ok]),
        "(c) both": set(trial_ids[site_ok & kin_ok]),
    }
    add(f"FC-valid liver sites           : {len(fc_sites):,}")
    add(f"Liver-graph protein nodes      : {len(liv_prot):,}")
    add(f"Trials rankable everywhere     : {len(keep_all):,}")
    for k, v in subsets.items():
        add(f"  {k:28s}: n = {len(v):,}")
    add("")

    # seed-averaged node2vec series
    def seed_avg(net: str, ids: set[int]) -> pd.Series | None:
        if len(n2v[net]) < 2:
            return None
        wide = None
        for s, d in n2v[net].items():
            dd = d[d[KEY].isin(ids)][[KEY, RANK_COL]].rename(columns={RANK_COL: f"s{s}"})
            wide = dd if wide is None else wide.merge(dd, on=KEY)
        cols = [c for c in wide.columns if c != KEY]
        return pd.Series(wide[cols].mean(axis=1).to_numpy(), index=wide[KEY].to_numpy())

    for label, ids in subsets.items():
        if len(ids) < 30:
            add(f"--- {label}: too few trials (n={len(ids)}), skipped ---")
            continue
        rows = []

        def add_test(da, db, la, lb):
            m = (da[da[KEY].isin(ids)][[KEY, RANK_COL]]
                 .merge(db[db[KEY].isin(ids)][[KEY, RANK_COL]], on=KEY,
                        suffixes=("_a", "_b")))
            a, b = m[f"{RANK_COL}_a"].to_numpy(), m[f"{RANK_COL}_b"].to_numpy()
            r = paired_test(a, b, la, lb)
            d = a - b
            r["B_wins"] = int((d > 0).sum()); r["A_wins"] = int((d < 0).sum())
            r["ties"] = int((d == 0).sum()); r["pct_ties"] = round(100.0 * (d == 0).mean(), 2)
            rows.append(r)

        add_test(gen_s, liv_s, "Generic (spec)", "Liver FC (spec)")
        for s, nd in sorted(nulls.items()):
            add_test(nd, liv_s, f"Null s{s} (spec)", "Liver FC (spec)")

        ga, la_ = seed_avg("generic", ids), seed_avg("liver_fc", ids)
        if ga is not None and la_ is not None:
            common = ga.index.intersection(la_.index)
            a, b = ga.loc[common].to_numpy(), la_.loc[common].to_numpy()
            r = paired_test(a, b, "Generic n2v seed-avg", "Liver FC n2v seed-avg")
            d = a - b
            r["B_wins"] = int((d > 0).sum()); r["A_wins"] = int((d < 0).sum())
            r["ties"] = int((d == 0).sum()); r["pct_ties"] = round(100.0 * (d == 0).mean(), 2)
            rows.append(r)

        t = pd.DataFrame(rows)
        t["p_holm"] = holm(t["p_value"].tolist())
        t["significant_holm_0.05"] = t["p_holm"] < 0.05
        cols = ["comparison", "n", "median_A", "median_B", "median_paired_diff",
                "rank_biserial", "B_wins", "A_wins", "ties", "pct_ties",
                "p_value", "p_holm", "significant_holm_0.05", "favours"]
        add(f"--- {label} (n = {len(ids):,}) ---")
        add(fmt(t[cols]))
        add("")
        t[cols].assign(subset=label).to_csv(
            OUT / f"subset_{re.sub(r'[^a-z0-9]+', '_', label.lower()).strip('_')}.csv",
            index=False)

    text = "\n".join(L)
    (OUT / "subset_analysis.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
