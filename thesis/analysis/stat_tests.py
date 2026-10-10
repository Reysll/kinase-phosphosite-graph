"""
G1 -- Significance testing for the four-network / two-embedding LOO comparison.

READ-ONLY. This script reads existing experiment outputs under
outputs_prediction/ and writes only into thesis/analysis/out/. It does not
modify any project code, data or results.

Comparison family (10 tests, Holm-corrected as ONE family)
----------------------------------------------------------
Block A -- context effect, within each embedding:
    Generic vs Control, Generic vs Cancer, Generic vs Liver FC   (x node2vec, spectral) = 6
Block B -- embedding effect, within each network:
    node2vec vs spectral                                          (x 4 networks)         = 4

Statistics per test
-------------------
  * paired Wilcoxon signed-rank on per-trial adjusted_held_out_rank
  * median paired difference (rank_A - rank_B; positive => B ranks better)
  * matched-pairs rank-biserial correlation (Kerby 2014)
  * n used per test

Analyses
--------
  1. trial-level, all rankable trials
  2. kinase-level sensitivity (mean adjusted rank per held-out kinase, then paired test)
  3. trial-level restricted to `informative` kinases (G15)

Unrankable trials
-----------------
The UNION of trials whose adjusted_held_out_rank is NaN in ANY of the eight
experiments is removed from EVERY comparison, so all tests use an identical
trial set.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.stat_tests

Seed: none required (all tests are deterministic given the inputs).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PRED = Path("outputs_prediction")
OUT = Path("thesis/analysis/out")

NETWORKS = {
    "Generic": "generic_multi_kinase",
    "Control": "control_multi_kinase",
    "Cancer": "cancer_multi_kinase",
    "Liver FC": "liver_multi_kinase",
}
EMBEDDINGS = {"node2vec": "", "spectral": "_spectral"}
RANK_COL = "adjusted_held_out_rank"
KEY = "trial_index"


# ----------------------------------------------------------------------------- loading
def results_path(network: str, embedding: str) -> Path:
    return PRED / f"{NETWORKS[network]}{EMBEDDINGS[embedding]}" / "results.csv.gz"


def git_hash_for(path: Path) -> str:
    """Last commit that touched *path* (provenance for the appendix)."""
    try:
        out = subprocess.run(
            ["git", "log", "-1", "--format=%H", "--", str(path)],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        return out or "UNTRACKED/UNCOMMITTED"
    except Exception:
        return "UNKNOWN"


def load_all() -> tuple[dict[tuple[str, str], pd.DataFrame], pd.DataFrame]:
    frames: dict[tuple[str, str], pd.DataFrame] = {}
    prov = []
    for emb in EMBEDDINGS:
        for net in NETWORKS:
            p = results_path(net, emb)
            df = pd.read_csv(p)[[KEY, "true_kinase_node_id", RANK_COL]]
            frames[(net, emb)] = df
            prov.append({
                "network": net, "embedding": emb, "path": str(p).replace("\\", "/"),
                "n_rows": len(df), "n_nan_rank": int(df[RANK_COL].isna().sum()),
                "git_commit": git_hash_for(p),
            })
    return frames, pd.DataFrame(prov)


def union_unrankable(frames: dict[tuple[str, str], pd.DataFrame]) -> set[int]:
    bad: set[int] = set()
    for df in frames.values():
        bad |= set(df.loc[df[RANK_COL].isna(), KEY].astype(int))
    return bad


# ----------------------------------------------------------------------------- stats
def rank_biserial(diff: np.ndarray) -> float:
    """Matched-pairs rank-biserial correlation (Kerby 2014).

    diff = rank_A - rank_B with zeros dropped (Wilcoxon convention).
    Positive result => B ranks better (lower) more often / more strongly.
    """
    d = diff[diff != 0]
    if d.size == 0:
        return float("nan")
    r = stats.rankdata(np.abs(d))
    w_pos = r[d > 0].sum()   # A worse  => favours B
    w_neg = r[d < 0].sum()   # A better => favours A
    total = w_pos + w_neg
    return float((w_pos - w_neg) / total) if total else float("nan")


def paired_test(a: np.ndarray, b: np.ndarray, label_a: str, label_b: str) -> dict:
    diff = a - b
    n_ties = int((diff == 0).sum())
    if np.all(diff == 0):
        stat, p = float("nan"), 1.0
    else:
        stat, p = stats.wilcoxon(a, b, zero_method="wilcox", alternative="two-sided")
    r_rb = rank_biserial(diff)
    med = float(np.median(diff))
    if np.isnan(r_rb) or med == 0:
        favours = "tie"
    else:
        favours = label_b if med > 0 else label_a
    return {
        "comparison": f"{label_a} vs {label_b}",
        "n": int(len(a)),
        "n_ties": n_ties,
        "median_A": float(np.median(a)),
        "median_B": float(np.median(b)),
        "median_paired_diff": med,
        "rank_biserial": r_rb,
        "wilcoxon_W": float(stat),
        "p_value": float(p),
        "favours": favours,
    }


def holm(pvals: list[float]) -> list[float]:
    """Holm-Bonferroni step-down adjusted p-values."""
    n = len(pvals)
    order = np.argsort(pvals)
    adj = np.empty(n, dtype=float)
    running = 0.0
    for i, idx in enumerate(order):
        val = (n - i) * pvals[idx]
        running = max(running, val)
        adj[idx] = min(1.0, running)
    return adj.tolist()


# ----------------------------------------------------------------------------- family
def build_family() -> list[tuple[str, str, str, str, str]]:
    """(block, label_A, emb_A, label_B, emb_B)"""
    fam = []
    for emb in ("node2vec", "spectral"):
        for other in ("Control", "Cancer", "Liver FC"):
            fam.append(("A: context", "Generic", emb, other, emb))
    for net in NETWORKS:
        fam.append(("B: embedding", net, "node2vec", net, "spectral"))
    return fam


def run_block(
    frames, keep_trials: set[int], kinase_filter: set[str] | None, level: str
) -> pd.DataFrame:
    """level: 'trial' or 'kinase'."""
    rows = []
    for block, net_a, emb_a, net_b, emb_b in build_family():
        da = frames[(net_a, emb_a)]
        db = frames[(net_b, emb_b)]
        da = da[da[KEY].isin(keep_trials)]
        db = db[db[KEY].isin(keep_trials)]
        if kinase_filter is not None:
            da = da[da["true_kinase_node_id"].isin(kinase_filter)]
            db = db[db["true_kinase_node_id"].isin(kinase_filter)]

        if level == "kinase":
            ga = da.groupby("true_kinase_node_id")[RANK_COL].mean()
            gb = db.groupby("true_kinase_node_id")[RANK_COL].mean()
            common = ga.index.intersection(gb.index)
            a, b = ga.loc[common].to_numpy(), gb.loc[common].to_numpy()
        else:
            m = da.merge(db, on=KEY, suffixes=("_a", "_b"))
            a, b = m[f"{RANK_COL}_a"].to_numpy(), m[f"{RANK_COL}_b"].to_numpy()

        lab_a = f"{net_a} ({emb_a})"
        lab_b = f"{net_b} ({emb_b})"
        res = paired_test(a, b, lab_a, lab_b)
        res["block"] = block
        rows.append(res)

    out = pd.DataFrame(rows)
    out["p_holm"] = holm(out["p_value"].tolist())
    out["significant_holm_0.05"] = out["p_holm"] < 0.05
    cols = ["block", "comparison", "n", "n_ties", "median_A", "median_B",
            "median_paired_diff", "rank_biserial", "wilcoxon_W", "p_value",
            "p_holm", "significant_holm_0.05", "favours"]
    return out[cols]


def informative_kinases(embedding: str) -> set[str]:
    p = PRED / "analysis" / f"rank_heatmap_{embedding}.csv"
    df = pd.read_csv(p)
    names = df.loc[df["informative"] == True, "kinase"].astype(str)  # noqa: E712
    return {f"PROTEIN:{n}" for n in names}


# ----------------------------------------------------------------------------- main
def fmt(df: pd.DataFrame) -> str:
    d = df.copy()
    for c in ("median_A", "median_B", "median_paired_diff", "wilcoxon_W"):
        d[c] = d[c].map(lambda v: f"{v:,.1f}")
    d["rank_biserial"] = d["rank_biserial"].map(lambda v: f"{v:+.3f}")
    for c in ("p_value", "p_holm"):
        d[c] = d[c].map(lambda v: f"{v:.3e}" if v < 1e-3 else f"{v:.4f}")
    return d.to_string(index=False)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames, prov = load_all()

    bad = union_unrankable(frames)
    all_trials = set(frames[("Generic", "node2vec")][KEY].astype(int))
    keep = all_trials - bad

    lines: list[str] = []
    add = lines.append
    add("G1 -- Significance testing, four-network / two-embedding LOO comparison")
    add("=" * 78)
    add("")
    add(f"Total trials in frozen set      : {len(all_trials):,}")
    add(f"Unrankable in >=1 experiment    : {len(bad):,}  (union across all 8 runs)")
    add(f"Trials used in EVERY comparison : {len(keep):,}")
    add("")
    add("Family: 10 paired tests, Holm-corrected together.")
    add("Convention: diff = rank_A - rank_B. Positive median_paired_diff and")
    add("positive rank_biserial both mean B achieves the BETTER (lower) rank.")
    add("")

    # 1. trial level
    t1 = run_block(frames, keep, None, "trial")
    add("1. TRIAL-LEVEL (primary)")
    add("-" * 78)
    add(fmt(t1))
    add("")

    # 2. kinase level
    t2 = run_block(frames, keep, None, "kinase")
    add("2. KINASE-LEVEL SENSITIVITY (mean adjusted rank per held-out kinase)")
    add("-" * 78)
    add(fmt(t2))
    add("")

    # 3. informative-restricted
    inf_n2v = informative_kinases("node2vec")
    inf_spec = informative_kinases("spectral")
    inf_both = inf_n2v & inf_spec
    add("3. INFORMATIVE-KINASE-RESTRICTED (G15), trial level")
    add("-" * 78)
    add(f"informative(node2vec) = {len(inf_n2v)}   informative(spectral) = {len(inf_spec)}"
        f"   intersection = {len(inf_both)}")
    add("Block A uses the tested embedding's own set; Block B uses the intersection.")
    add("")
    rows = []
    for block, na, ea, nb, eb in build_family():
        filt = inf_both if block.startswith("B") else informative_kinases(ea)
        sub = run_block(frames, keep, filt, "trial")
        rows.append(sub[sub["comparison"] == f"{na} ({ea}) vs {nb} ({eb})"].iloc[0])
    t3 = pd.DataFrame(rows).reset_index(drop=True)
    t3["p_holm"] = holm(t3["p_value"].tolist())
    t3["significant_holm_0.05"] = t3["p_holm"] < 0.05
    add(fmt(t3))
    add("")

    add("PROVENANCE (git commit of each results file)")
    add("-" * 78)
    add(prov.to_string(index=False))

    text = "\n".join(lines)
    (OUT / "stat_tests.txt").write_text(text, encoding="utf-8")
    t1.to_csv(OUT / "stat_tests_trial_level.csv", index=False)
    t2.to_csv(OUT / "stat_tests_kinase_level.csv", index=False)
    t3.to_csv(OUT / "stat_tests_informative.csv", index=False)
    prov.to_csv(OUT / "stat_tests_provenance.csv", index=False)
    print(text)


if __name__ == "__main__":
    main()
