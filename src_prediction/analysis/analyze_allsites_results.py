"""
Performance summary for the cluster all-sites spectral LOO runs.

Reads from: cluster_results/{network}_spectral_dim64_allsites/
Writes to:  outputs_prediction/analysis_allsites/

Does not touch any existing outputs_prediction/analysis/ files.
Produces two sections:
  1. All-sites (14,596 trials)
  2. Multi-kinase filtered (n_true_kinases_for_site >= 2, ~6,909 trials)
     for direct comparison with Phase 2.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

CLUSTER_RESULTS = Path("cluster_results")
OUT_DIR = Path("outputs_prediction/analysis_allsites")
OUT_DIR.mkdir(parents=True, exist_ok=True)

EXPERIMENTS = {
    "Generic": "generic_spectral_dim64_allsites",
    "Control": "control_spectral_dim64_allsites",
    "Cancer": "cancer_spectral_dim64_allsites",
    "Liver FC": "liver_fc_spectral_dim64_allsites",
}

# Phase 2 multi-kinase spectral results (6,909 trials, dim=32) for comparison
PHASE2_MEDIAN = {"Generic": 79, "Control": 74, "Cancer": 74, "Liver FC": 74}
PHASE2_TOP20 = {"Generic": 10.6, "Control": 22.1, "Cancer": 22.1, "Liver FC": 22.0}


def load_metrics(exp_dir: Path) -> dict:
    path = exp_dir / "metrics.csv.gz"
    df = pd.read_csv(path).set_index("metric")["value"]
    return df.to_dict()


def compute_metrics_from_results(df: pd.DataFrame) -> dict:
    """Recompute key metrics from a results DataFrame (supports filtered subsets)."""
    ranked = df.dropna(subset=["adjusted_held_out_rank"])
    n = len(ranked)
    adj = ranked["adjusted_held_out_rank"]
    best = ranked["best_true_kinase_rank"].dropna()
    return {
        "adjusted_held_out_median_rank": adj.median(),
        "adjusted_held_out_mean_rank": adj.mean(),
        "adjusted_held_out_top_1_accuracy": (adj <= 1).mean(),
        "adjusted_held_out_top_5_accuracy": (adj <= 5).mean(),
        "adjusted_held_out_top_10_accuracy": (adj <= 10).mean(),
        "adjusted_held_out_top_20_accuracy": (adj <= 20).mean(),
        "best_true_median_rank": best.median(),
        "best_true_top_20_accuracy": (best <= 20).mean(),
        "adjusted_held_out_n_ranked": float(n),
    }


def load_results(exp_dir: Path) -> pd.DataFrame:
    return pd.read_csv(exp_dir / "results.csv.gz")


def fmt_pct(v: float) -> str:
    return f"{v * 100:.1f}%"


def main() -> None:
    rows = {}
    for label, dirname in EXPERIMENTS.items():
        exp_dir = CLUSTER_RESULTS / dirname
        if not exp_dir.exists():
            print(f"WARNING: missing {exp_dir}")
            continue
        rows[label] = load_metrics(exp_dir)

    networks = list(rows.keys())

    lines = []
    lines.append("All-Sites Spectral LOO Results (Cluster Run)")
    lines.append("=" * 70)
    lines.append(f"Trials:     14,596  (all sites: single-kinase + multi-kinase)")
    lines.append(f"Embedding:  Spectral, dim=64")
    lines.append(f"Candidates: 420 kinases")
    lines.append("")

    # Header
    col_w = 14
    lines.append(f"{'Metric':<32}" + "".join(f"{n:>{col_w}}" for n in networks))
    lines.append("-" * (32 + col_w * len(networks)))

    def row(label, key, fmt=lambda x: f"{x:.1f}"):
        vals = "".join(f"{fmt(rows[n][key]):>{col_w}}" for n in networks)
        lines.append(f"{label:<32}{vals}")

    row("Adj median rank (lower=better)", "adjusted_held_out_median_rank")
    row("Adj mean rank", "adjusted_held_out_mean_rank")
    row("Adj top-1 %", "adjusted_held_out_top_1_accuracy", fmt_pct)
    row("Adj top-5 %", "adjusted_held_out_top_5_accuracy", fmt_pct)
    row("Adj top-10 %", "adjusted_held_out_top_10_accuracy", fmt_pct)
    row("Adj top-20 %", "adjusted_held_out_top_20_accuracy", fmt_pct)
    row("Best-true median rank", "best_true_median_rank")
    row("Best-true top-20 %", "best_true_top_20_accuracy", fmt_pct)
    row("N trials ranked", "adjusted_held_out_n_ranked", fmt=lambda x: f"{int(x):,}")

    lines.append("")
    lines.append("Comparison vs Phase 2 multi-kinase spectral (6,909 trials, dim=32)")
    lines.append("-" * (32 + col_w * len(networks)))
    lines.append(f"{'':32}" + "".join(f"{n:>{col_w}}" for n in networks))

    p2_med = "".join(f"{PHASE2_MEDIAN[n]:>{col_w}}" for n in networks)
    lines.append(f"{'Phase 2 adj median rank':<32}{p2_med}")

    allsites_med = "".join(
        f"{int(rows[n]['adjusted_held_out_median_rank']):>{col_w}}" for n in networks
    )
    lines.append(f"{'All-sites adj median rank':<32}{allsites_med}")

    delta = "".join(
        f"{int(rows[n]['adjusted_held_out_median_rank']) - PHASE2_MEDIAN[n]:>+{col_w}}"
        for n in networks
    )
    lines.append(f"{'Delta (all-sites minus Phase 2)':<32}{delta}")

    lines.append("")
    p2_top20 = "".join(f"{PHASE2_TOP20[n]:>{col_w - 1}.1f}%" for n in networks)
    lines.append(f"{'Phase 2 adj top-20 %':<32}{p2_top20}")

    allsites_top20 = "".join(
        f"{rows[n]['adjusted_held_out_top_20_accuracy'] * 100:>{col_w - 1}.1f}%" for n in networks
    )
    lines.append(f"{'All-sites adj top-20 %':<32}{allsites_top20}")

    lines.append("")
    lines.append("Notes")
    lines.append("-----")
    lines.append("All-sites includes 7,687 single-kinase sites (76% of trials). Single-kinase")
    lines.append("LOO is less informative: no within-site negative contrast (only one true")
    lines.append("kinase, so adjusted rank = held-out rank). Multi-kinase median rank remains")
    lines.append("the cleaner metric for cross-network comparison.")
    lines.append("dim=64 here vs dim=32 in Phase 2; higher dimensions may contribute a small")
    lines.append("fraction of the performance difference.")

    # ------------------------------------------------------------------
    # Section 2: multi-kinase filtered (n_true_kinases_for_site >= 2)
    # ------------------------------------------------------------------
    print("Loading results for multi-kinase filter (this may take a moment)...")
    multi_rows = {}
    for label, dirname in EXPERIMENTS.items():
        exp_dir = CLUSTER_RESULTS / dirname
        if not exp_dir.exists():
            continue
        df = load_results(exp_dir)
        multi = df[df["n_true_kinases_for_site"] >= 2]
        multi_rows[label] = compute_metrics_from_results(multi)

    n_multi = int(next(iter(multi_rows.values()))["adjusted_held_out_n_ranked"])

    lines.append("")
    lines.append("")
    lines.append("Multi-Kinase Filtered Results (n_true_kinases >= 2)")
    lines.append("=" * 70)
    lines.append(f"Trials:     ~{n_multi:,}  (sites with 2+ known kinases only)")
    lines.append(f"Embedding:  Spectral, dim=64  |  graph: all-sites embedding")
    lines.append(f"Candidates: 420 kinases")
    lines.append("")
    lines.append(f"{'Metric':<32}" + "".join(f"{n:>{col_w}}" for n in networks))
    lines.append("-" * (32 + col_w * len(networks)))

    def mrow(label, key, fmt=lambda x: f"{x:.1f}"):
        vals = "".join(f"{fmt(multi_rows[n][key]):>{col_w}}" for n in networks)
        lines.append(f"{label:<32}{vals}")

    mrow("Adj median rank (lower=better)", "adjusted_held_out_median_rank")
    mrow("Adj mean rank", "adjusted_held_out_mean_rank")
    mrow("Adj top-1 %", "adjusted_held_out_top_1_accuracy", fmt_pct)
    mrow("Adj top-5 %", "adjusted_held_out_top_5_accuracy", fmt_pct)
    mrow("Adj top-10 %", "adjusted_held_out_top_10_accuracy", fmt_pct)
    mrow("Adj top-20 %", "adjusted_held_out_top_20_accuracy", fmt_pct)
    mrow("Best-true median rank", "best_true_median_rank")
    mrow("Best-true top-20 %", "best_true_top_20_accuracy", fmt_pct)
    mrow("N trials ranked", "adjusted_held_out_n_ranked", fmt=lambda x: f"{int(x):,}")

    lines.append("")
    lines.append("Direct comparison vs Phase 2 (same trial type, same filter)")
    lines.append("-" * (32 + col_w * len(networks)))
    lines.append(f"{'':32}" + "".join(f"{n:>{col_w}}" for n in networks))

    p2_med = "".join(f"{PHASE2_MEDIAN[n]:>{col_w}}" for n in networks)
    lines.append(f"{'Phase 2 adj median rank':<32}{p2_med}")

    multi_med = "".join(
        f"{int(multi_rows[n]['adjusted_held_out_median_rank']):>{col_w}}" for n in networks
    )
    lines.append(f"{'Cluster multi-filt median rank':<32}{multi_med}")

    delta2 = "".join(
        f"{int(multi_rows[n]['adjusted_held_out_median_rank']) - PHASE2_MEDIAN[n]:>+{col_w}}"
        for n in networks
    )
    lines.append(f"{'Delta (cluster minus Phase 2)':<32}{delta2}")

    lines.append("")
    p2_top20 = "".join(f"{PHASE2_TOP20[n]:>{col_w - 1}.1f}%" for n in networks)
    lines.append(f"{'Phase 2 adj top-20 %':<32}{p2_top20}")

    multi_top20 = "".join(
        f"{multi_rows[n]['adjusted_held_out_top_20_accuracy'] * 100:>{col_w - 1}.1f}%"
        for n in networks
    )
    lines.append(f"{'Cluster multi-filt top-20 %':<32}{multi_top20}")

    lines.append("")
    lines.append("Note: cluster embedding trained on all-sites graph (dim=64);")
    lines.append("Phase 2 embedding trained on multi-kinase-only graph (dim=32).")
    lines.append("Remaining delta reflects dim difference and embedding scope.")

    summary = "\n".join(lines)
    out_path = OUT_DIR / "allsites_spectral_summary.txt"
    out_path.write_text(summary)
    print(summary)
    print(f"\nWrote: {out_path}")


if __name__ == "__main__":
    main()
