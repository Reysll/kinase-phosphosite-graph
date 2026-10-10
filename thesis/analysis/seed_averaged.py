"""
Item 2 -- seed-averaged node2vec analysis (read-only).

node2vec results depend on the random seed. Rather than reporting one seed, this
averages each trial's adjusted_held_out_rank across seeds {42, 7, 13} within a
network, then runs the Generic vs Liver FC contrast on those averaged ranks.
Averaging suppresses per-seed sampling noise, so the surviving difference is the
part of the context effect that is stable across embedding draws.

Also reports, per network:
    * per-seed median adjusted rank
    * cross-seed spread (max - min of the per-seed medians) and SD
    * per-trial cross-seed SD, summarised
    * top-1 predicted kinase monopoly per seed (which hub wins, and how often)

Inputs (must exist):
    outputs_prediction/{generic,liver}_multi_kinase/results.csv.gz       seed 42
    thesis/analysis/out/seed_variance/seedvar_{net}_node2vec_seed{7,13}_results.csv.gz

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.seed_averaged
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from thesis.analysis.stat_tests import KEY, RANK_COL, holm, paired_test
from thesis.analysis.verify_tests import fmt

PRED = Path("outputs_prediction")
SEED_DIR = Path("thesis/analysis/out/seed_variance")
OUT = Path("thesis/analysis/out")
SEEDS = [42, 7, 13]

BASE42 = {
    "generic": PRED / "generic_multi_kinase" / "results.csv.gz",
    "liver_fc": PRED / "liver_multi_kinase" / "results.csv.gz",
}


def path_for(net: str, seed: int) -> Path:
    if seed == 42:
        return BASE42[net]
    return SEED_DIR / f"seedvar_{net}_node2vec_seed{seed}_results.csv.gz"


def load(net: str, seed: int) -> pd.DataFrame:
    df = pd.read_csv(path_for(net, seed))
    keep = [KEY, "true_kinase_node_id", RANK_COL]
    if "top1_predicted_kinase" in df.columns:
        keep.append("top1_predicted_kinase")
    return df[keep]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L: list[str] = ["ITEM 2 -- SEED-AVERAGED node2vec ANALYSIS", "=" * 78, ""]
    add = L.append

    frames: dict[tuple[str, int], pd.DataFrame] = {}
    missing = []
    for net in BASE42:
        for s in SEEDS:
            p = path_for(net, s)
            if not p.exists():
                missing.append(str(p))
            else:
                frames[(net, s)] = load(net, s)
    if missing:
        add("MISSING INPUTS:")
        for m in missing:
            add("  " + m)
        print("\n".join(L))
        return

    # trials unrankable in ANY run are dropped everywhere
    bad: set[int] = set()
    for df in frames.values():
        bad |= set(df.loc[df[RANK_COL].isna(), KEY].astype(int))
    keep = set(frames[("generic", 42)][KEY].astype(int)) - bad
    add(f"Seeds averaged         : {SEEDS}")
    add(f"Trials dropped (NaN in any run): {len(bad)}")
    add(f"Trials used            : {len(keep):,}")
    add("")

    # ---------------------------------------------------------- per-seed medians
    add("PER-SEED MEDIAN ADJUSTED RANK, AND CROSS-SEED SPREAD")
    add("-" * 78)
    rows = []
    for net in BASE42:
        meds = {}
        for s in SEEDS:
            d = frames[(net, s)]
            meds[f"seed{s}"] = float(d[d[KEY].isin(keep)][RANK_COL].median())
        vals = np.array(list(meds.values()), dtype=float)
        rows.append({"network": net, **meds,
                     "spread_max_minus_min": float(vals.max() - vals.min()),
                     "sd_of_medians": float(vals.std(ddof=1))})
    med_t = pd.DataFrame(rows)
    add(med_t.to_string(index=False))
    med_t.to_csv(OUT / "seed_averaged_medians.csv", index=False)
    add("")

    # ---------------------------------------------------------- per-trial spread
    add("PER-TRIAL CROSS-SEED VARIABILITY (SD of a trial's rank across the 3 seeds)")
    add("-" * 78)
    avg: dict[str, pd.Series] = {}
    for net in BASE42:
        wide = None
        for s in SEEDS:
            d = frames[(net, s)]
            d = d[d[KEY].isin(keep)][[KEY, RANK_COL]].rename(columns={RANK_COL: f"s{s}"})
            wide = d if wide is None else wide.merge(d, on=KEY)
        wide = wide.sort_values(KEY).reset_index(drop=True)
        cols = [f"s{s}" for s in SEEDS]
        sd = wide[cols].std(axis=1, ddof=1)
        avg[net] = pd.Series(wide[cols].mean(axis=1).to_numpy(), index=wide[KEY].to_numpy())
        add(f"  {net:9s}  mean SD={sd.mean():6.1f}   median SD={sd.median():6.1f}   "
            f"p90 SD={sd.quantile(0.9):6.1f}   max SD={sd.max():6.1f}")
    add("")

    # ---------------------------------------------------------- top-1 monopoly
    add("TOP-1 PREDICTED KINASE PER SEED (hub identity stability)")
    add("-" * 78)
    t1rows = []
    for net in BASE42:
        for s in SEEDS:
            d = frames[(net, s)]
            d = d[d[KEY].isin(keep)]
            if "top1_predicted_kinase" not in d.columns:
                continue
            vc = d["top1_predicted_kinase"].value_counts()
            t1rows.append({
                "network": net, "seed": s,
                "top1_kinase": str(vc.index[0]).split(":")[-1],
                "pct_of_trials": round(100.0 * vc.iloc[0] / len(d), 1),
                "n_distinct_top1": int(d["top1_predicted_kinase"].nunique()),
            })
    if t1rows:
        t1 = pd.DataFrame(t1rows)
        add(t1.to_string(index=False))
        t1.to_csv(OUT / "seed_averaged_top1.csv", index=False)
    add("")

    # ---------------------------------------------------------- the contrast
    add("GENERIC vs LIVER FC ON SEED-AVERAGED RANKS")
    add("-" * 78)
    common = avg["generic"].index.intersection(avg["liver_fc"].index)
    a = avg["generic"].loc[common].to_numpy()
    b = avg["liver_fc"].loc[common].to_numpy()
    res = paired_test(a, b, "Generic (seed-avg)", "Liver FC (seed-avg)")
    d = a - b
    res["B_wins"] = int((d > 0).sum())
    res["A_wins"] = int((d < 0).sum())
    res["ties"] = int((d == 0).sum())
    res["pct_ties"] = round(100.0 * (d == 0).mean(), 2)

    # single-seed references, same trials
    ref = []
    for s in SEEDS:
        ga = frames[("generic", s)]
        lb = frames[("liver_fc", s)]
        m = (ga[ga[KEY].isin(keep)][[KEY, RANK_COL]]
             .merge(lb[lb[KEY].isin(keep)][[KEY, RANK_COL]], on=KEY,
                    suffixes=("_g", "_l")))
        r = paired_test(m[f"{RANK_COL}_g"].to_numpy(), m[f"{RANK_COL}_l"].to_numpy(),
                        f"Generic seed{s}", f"Liver FC seed{s}")
        dd = m[f"{RANK_COL}_g"].to_numpy() - m[f"{RANK_COL}_l"].to_numpy()
        r["B_wins"] = int((dd > 0).sum()); r["A_wins"] = int((dd < 0).sum())
        r["ties"] = int((dd == 0).sum()); r["pct_ties"] = round(100.0 * (dd == 0).mean(), 2)
        ref.append(r)

    allt = pd.DataFrame([res] + ref)
    allt["p_holm"] = holm(allt["p_value"].tolist())
    allt["significant_holm_0.05"] = allt["p_holm"] < 0.05
    cols = ["comparison", "n", "median_A", "median_B", "median_paired_diff",
            "rank_biserial", "B_wins", "A_wins", "ties", "pct_ties",
            "p_value", "p_holm", "significant_holm_0.05", "favours"]
    add(fmt(allt[cols]))
    allt[cols].to_csv(OUT / "seed_averaged_contrast.csv", index=False)

    text = "\n".join(L)
    (OUT / "seed_averaged.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
