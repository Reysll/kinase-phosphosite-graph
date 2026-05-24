"""
Four supervisor-facing analyses over the multi-kinase LOO results.

Outputs written to outputs_prediction/analysis/:
  performance_summary.txt        — all four experiments side-by-side
  top1_frequency_generic.csv     — top-1 predicted kinase counts (generic model)
  top1_frequency_liver.csv       — top-1 predicted kinase counts (liver model)
  top1_frequency_comparison.csv  — side-by-side with rank and count delta
  top1_frequency_summary.txt     — human-readable top-20 for each model
  ttk_deepdive.txt               — TTK-specific accuracy and rank metrics
  ranking_shifts.txt             — top-20 biggest improvements and worsenings in liver
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src_prediction.config import PRED_OUTPUTS_DIR
from src_prediction.io_utils import write_text


ANALYSIS_DIR = PRED_OUTPUTS_DIR / "analysis"

EXPERIMENT_DIRS = {
    "generic_debug": PRED_OUTPUTS_DIR / "generic_trained_model_debug",
    "liver_debug": PRED_OUTPUTS_DIR / "liver_trained_model_fc_corr_debug",
    "generic_multi": PRED_OUTPUTS_DIR / "generic_multi_kinase",
    "liver_multi": PRED_OUTPUTS_DIR / "liver_multi_kinase",
}

EXPERIMENT_LABELS = {
    "generic_debug": "Generic model     (50-trial single-kinase debug)",
    "liver_debug": "Liver model       (50-trial single-kinase debug)",
    "generic_multi": "Generic model     (6,909-trial multi-kinase)",
    "liver_multi": "Liver model       (6,909-trial multi-kinase)",
}

TTK_NODE_ID = "PROTEIN:TTK"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _load_results(exp_key: str) -> pd.DataFrame | None:
    path = EXPERIMENT_DIRS[exp_key] / "results.csv.gz"
    if not path.exists():
        return None
    return pd.read_csv(path)


def _clean_node_id(node_id: str) -> str:
    """PROTEIN:GENE -> GENE,  SITE:GENE_pSXXX -> GENE_pSXXX"""
    if isinstance(node_id, str):
        return node_id.split(":", 1)[-1]
    return str(node_id)


def _top_k_str(series: pd.Series, k: int, n: int) -> str:
    valid = series.dropna()
    count = int((valid <= k).sum())
    pct = count / n if n else 0.0
    return f"{pct:.1%}  ({count}/{n})"


def _rank_stats(series: pd.Series, n: int) -> tuple[str, str]:
    valid = series.dropna()
    if valid.empty:
        return "N/A", "N/A"
    return f"{valid.mean():.1f}", f"{valid.median():.0f}"


# ---------------------------------------------------------------------------
# Analysis 1: combined performance summary
# ---------------------------------------------------------------------------


def analysis_performance_summary() -> str:
    col_w = 52

    lines = [
        "Combined Performance Summary — LOO Evaluation",
        "=" * 90,
        "",
        f"{'Experiment':<{col_w}}  {'n_trials':>8}  {'Top-1':>7}  {'Top-5':>7}  {'Top-10':>7}  {'Top-20':>7}  {'Mean rank':>10}  {'Median rank':>11}",
        "-" * 90,
        "Held-out kinase rank (the specific removed edge)",
        "-" * 90,
    ]

    for key, label in EXPERIMENT_LABELS.items():
        df = _load_results(key)
        if df is None:
            lines.append(f"  {label}  —  results not found")
            continue
        n = len(df)
        col = "held_out_kinase_rank"
        mean_r, med_r = _rank_stats(df[col], n)
        lines.append(
            f"  {label:<{col_w - 2}}"
            f"  {n:>8,}"
            f"  {_top_k_str(df[col], 1, n):>7}"
            f"  {_top_k_str(df[col], 5, n):>7}"
            f"  {_top_k_str(df[col], 10, n):>7}"
            f"  {_top_k_str(df[col], 20, n):>7}"
            f"  {mean_r:>10}"
            f"  {med_r:>11}"
        )

    lines += [
        "",
        "-" * 90,
        "Best-true kinase rank (best-ranked among all known kinases for the site)",
        "-" * 90,
    ]

    for key, label in EXPERIMENT_LABELS.items():
        df = _load_results(key)
        if df is None:
            continue
        n = len(df)
        col = "best_true_kinase_rank"
        if col not in df.columns:
            lines.append(f"  {label:<{col_w - 2}}  — column not present")
            continue
        mean_r, med_r = _rank_stats(df[col], n)
        lines.append(
            f"  {label:<{col_w - 2}}"
            f"  {n:>8,}"
            f"  {_top_k_str(df[col], 1, n):>7}"
            f"  {_top_k_str(df[col], 5, n):>7}"
            f"  {_top_k_str(df[col], 10, n):>7}"
            f"  {_top_k_str(df[col], 20, n):>7}"
            f"  {mean_r:>10}"
            f"  {med_r:>11}"
        )

    lines += [
        "",
        "Note: candidate kinase pool = 420. Random top-1 baseline ~0.24%, top-20 ~4.8%.",
        "Note: multi-kinase generic outperforms liver on held-out rank because multi-kinase",
        "      sites are dominated by well-studied hub kinases (AKT1, GSK3B, TP53) that are",
        "      richly connected in the generic graph. Liver wins for sparse, less-studied sites.",
    ]

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 2: top-1 kinase frequency
# ---------------------------------------------------------------------------


def analysis_top1_frequency() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, str]:
    generic_df = _load_results("generic_multi")
    liver_df = _load_results("liver_multi")

    def freq_table(df: pd.DataFrame, model: str) -> pd.DataFrame:
        counts = (
            df["top1_predicted_kinase"].dropna().apply(_clean_node_id).value_counts().reset_index()
        )
        counts.columns = ["kinase", f"top1_count_{model}"]
        counts[f"top1_pct_{model}"] = counts[f"top1_count_{model}"] / len(df)
        counts[f"rank_{model}"] = range(1, len(counts) + 1)
        return counts

    gen_freq = freq_table(generic_df, "generic")
    liv_freq = freq_table(liver_df, "liver")

    comparison = gen_freq.merge(liv_freq, on="kinase", how="outer").fillna(0)
    comparison["top1_count_generic"] = comparison["top1_count_generic"].astype(int)
    comparison["top1_count_liver"] = comparison["top1_count_liver"].astype(int)
    comparison["count_delta_liver_minus_generic"] = (
        comparison["top1_count_liver"] - comparison["top1_count_generic"]
    )
    comparison = comparison.sort_values("top1_count_liver", ascending=False).reset_index(drop=True)

    n_gen = len(generic_df)
    n_liv = len(liver_df)
    top20_gen = gen_freq.head(20)
    top20_liv = liv_freq.head(20)

    lines = [
        "Top-1 Predicted Kinase Frequency — Multi-Kinase LOO (6,909 trials per model)",
        "=" * 80,
        "",
        f"{'Rank':<6}  {'Generic model top-1 kinase':<28}  {'Count':>6}  {'%':>6}",
        "-" * 55,
    ]
    for _, row in top20_gen.iterrows():
        lines.append(
            f"  {int(row['rank_generic']):<4}  {row['kinase']:<28}  {int(row['top1_count_generic']):>6}  {row['top1_pct_generic']:>5.1%}"
        )

    lines += [
        "",
        f"{'Rank':<6}  {'Liver model top-1 kinase':<28}  {'Count':>6}  {'%':>6}",
        "-" * 55,
    ]
    for _, row in top20_liv.iterrows():
        lines.append(
            f"  {int(row['rank_liver']):<4}  {row['kinase']:<28}  {int(row['top1_count_liver']):>6}  {row['top1_pct_liver']:>5.1%}"
        )

    # Kinases that gained most in liver
    gained = comparison[comparison["count_delta_liver_minus_generic"] > 0].nlargest(
        10, "count_delta_liver_minus_generic"
    )
    lost = comparison[comparison["count_delta_liver_minus_generic"] < 0].nsmallest(
        10, "count_delta_liver_minus_generic"
    )

    lines += [
        "",
        "Top-10 kinases gaining most top-1 predictions in liver vs generic:",
        f"  {'Kinase':<28}  {'Generic':>8}  {'Liver':>8}  {'Delta':>8}",
        "  " + "-" * 52,
    ]
    for _, row in gained.iterrows():
        lines.append(
            f"  {row['kinase']:<28}  {int(row['top1_count_generic']):>8}  {int(row['top1_count_liver']):>8}  +{int(row['count_delta_liver_minus_generic']):>7}"
        )

    lines += [
        "",
        "Top-10 kinases losing top-1 predictions in liver vs generic:",
        f"  {'Kinase':<28}  {'Generic':>8}  {'Liver':>8}  {'Delta':>8}",
        "  " + "-" * 52,
    ]
    for _, row in lost.iterrows():
        lines.append(
            f"  {row['kinase']:<28}  {int(row['top1_count_generic']):>8}  {int(row['top1_count_liver']):>8}  {int(row['count_delta_liver_minus_generic']):>8}"
        )

    # TTK specifically
    ttk_gen = int(gen_freq.loc[gen_freq["kinase"] == "TTK", "top1_count_generic"].sum())
    ttk_liv = int(liv_freq.loc[liv_freq["kinase"] == "TTK", "top1_count_liver"].sum())
    lines += [
        "",
        f"TTK top-1 predictions: generic = {ttk_gen}  |  liver = {ttk_liv}  (out of 6,909 trials each)",
    ]

    return gen_freq, liv_freq, comparison, "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 3: TTK deep-dive
# ---------------------------------------------------------------------------


def analysis_ttk_deepdive() -> str:
    lines = [
        "TTK Deep-Dive — LOO Evaluation across all TTK substrates",
        "=" * 70,
        "",
        "TTK is a liver cancer-relevant kinase with 68 substrates in the graph.",
        "",
    ]

    for key, label in EXPERIMENT_LABELS.items():
        df = _load_results(key)
        if df is None:
            lines.append(f"{label}  — results not found\n")
            continue

        ttk_rows = df[df["true_kinase_node_id"] == TTK_NODE_ID].copy()
        n = len(ttk_rows)

        lines.append(f"{label}")
        lines.append("-" * len(label))

        if n == 0:
            lines.append("  No TTK trials found in this experiment.\n")
            continue

        ho = ttk_rows["held_out_kinase_rank"]
        bt = (
            ttk_rows["best_true_kinase_rank"]
            if "best_true_kinase_rank" in ttk_rows.columns
            else pd.Series(dtype=float)
        )

        lines.append(f"  TTK trials (sites where TTK is the held-out kinase):  {n}")
        lines.append(
            f"  Held-out rank:    top-1 {_top_k_str(ho, 1, n)}  top-5 {_top_k_str(ho, 5, n)}  top-10 {_top_k_str(ho, 10, n)}  top-20 {_top_k_str(ho, 20, n)}"
        )
        lines.append(
            f"                    mean {ho.mean():.1f}  median {ho.median():.0f}  min {ho.min():.0f}  max {ho.max():.0f}"
        )

        if not bt.empty and bt.notna().any():
            lines.append(
                f"  Best-true rank:   top-1 {_top_k_str(bt, 1, n)}  top-5 {_top_k_str(bt, 5, n)}  top-10 {_top_k_str(bt, 10, n)}  top-20 {_top_k_str(bt, 20, n)}"
            )
            lines.append(f"                    mean {bt.mean():.1f}  median {bt.median():.0f}")

        # Top-1 predictions when TTK is the true kinase
        top1_counts = (
            ttk_rows["top1_predicted_kinase"].dropna().apply(_clean_node_id).value_counts().head(10)
        )
        lines.append(f"  Top-10 most frequent top-1 predictions (when TTK is true kinase):")
        for kinase, count in top1_counts.items():
            marker = " <-- TTK predicted correctly" if kinase == "TTK" else ""
            lines.append(f"    {kinase:<28}  {count:>4} times{marker}")

        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 4: ranking shift case studies
# ---------------------------------------------------------------------------


def analysis_ranking_shifts() -> str:
    comparison_path = PRED_OUTPUTS_DIR / "comparison_multi_kinase_generic_vs_liver.csv.gz"
    if not comparison_path.exists():
        return "Run compare_multi_kinase.py first to generate the comparison file."

    merged = pd.read_csv(comparison_path)

    def fmt_row(row: pd.Series) -> str:
        site = _clean_node_id(str(row["site_node_id"]))
        kinase = _clean_node_id(str(row["true_kinase_node_id"]))
        g_rank = (
            f"{int(row['generic_held_out_rank'])}"
            if pd.notna(row["generic_held_out_rank"])
            else "N/A"
        )
        l_rank = (
            f"{int(row['liver_held_out_rank'])}" if pd.notna(row["liver_held_out_rank"]) else "N/A"
        )
        delta = f"{row['held_out_rank_delta']:+.0f}"
        g_top1 = _clean_node_id(str(row.get("generic_top1", "?")))
        l_top1 = _clean_node_id(str(row.get("liver_top1", "?")))
        return (
            f"  {site:<28}  true: {kinase:<16}  "
            f"generic rank: {g_rank:>4}  liver rank: {l_rank:>4}  delta: {delta:>6}  "
            f"top-1: {g_top1} -> {l_top1}"
        )

    improved = merged[merged["held_out_rank_delta"] > 0].nlargest(20, "held_out_rank_delta")
    worsened = merged[merged["held_out_rank_delta"] < 0].nsmallest(20, "held_out_rank_delta")

    lines = [
        "Ranking Shift Case Studies — Multi-Kinase LOO (liver vs generic)",
        "=" * 90,
        "",
        "Interpretation: delta = generic_rank - liver_rank",
        "  Positive delta -> liver improved ranking (lower rank number = better)",
        "  Negative delta -> liver worsened ranking",
        "",
        f"Top-20 trials with BIGGEST IMPROVEMENT in liver model (n_total_improved = {int(merged['held_out_improved_in_liver'].sum()):,}):",
        "-" * 90,
    ]
    for _, row in improved.iterrows():
        lines.append(fmt_row(row))

    lines += [
        "",
        f"Top-20 trials with BIGGEST WORSENING in liver model (n_total_worsened = {int((~merged['held_out_improved_in_liver']).sum()):,}):",
        "-" * 90,
    ]
    for _, row in worsened.iterrows():
        lines.append(fmt_row(row))

    # Summarise by true kinase: which kinases improve/worsen most in liver?
    kinase_summary = (
        merged.groupby("true_kinase_node_id")
        .agg(
            n_trials=("trial_index", "count"),
            mean_delta=("held_out_rank_delta", "mean"),
            n_improved=("held_out_improved_in_liver", "sum"),
        )
        .reset_index()
    )
    kinase_summary["kinase"] = kinase_summary["true_kinase_node_id"].apply(_clean_node_id)
    kinase_summary["pct_improved"] = kinase_summary["n_improved"] / kinase_summary["n_trials"]
    kinase_summary = kinase_summary[kinase_summary["n_trials"] >= 5]  # only kinases with >=5 trials

    top_improved_kinases = kinase_summary.nlargest(15, "mean_delta")
    top_worsened_kinases = kinase_summary.nsmallest(15, "mean_delta")

    lines += [
        "",
        "Kinases with most consistent IMPROVEMENT in liver (min 5 trials, sorted by mean delta):",
        f"  {'Kinase':<20}  {'n_trials':>8}  {'mean_delta':>10}  {'n_improved':>10}  {'pct_improved':>12}",
        "  " + "-" * 64,
    ]
    for _, row in top_improved_kinases.iterrows():
        lines.append(
            f"  {row['kinase']:<20}  {int(row['n_trials']):>8}  {row['mean_delta']:>+10.1f}  {int(row['n_improved']):>10}  {row['pct_improved']:>11.1%}"
        )

    lines += [
        "",
        "Kinases with most consistent WORSENING in liver (min 5 trials):",
        f"  {'Kinase':<20}  {'n_trials':>8}  {'mean_delta':>10}  {'n_improved':>10}  {'pct_improved':>12}",
        "  " + "-" * 64,
    ]
    for _, row in top_worsened_kinases.iterrows():
        lines.append(
            f"  {row['kinase']:<20}  {int(row['n_trials']):>8}  {row['mean_delta']:>+10.1f}  {int(row['n_improved']):>10}  {row['pct_improved']:>11.1%}"
        )

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Analysis 1: Performance summary ===")
    perf = analysis_performance_summary()
    write_text(perf, ANALYSIS_DIR / "performance_summary.txt")
    print(perf)
    print()

    print("=== Analysis 2: Top-1 kinase frequency ===")
    gen_freq, liv_freq, comparison, freq_text = analysis_top1_frequency()
    gen_freq.to_csv(ANALYSIS_DIR / "top1_frequency_generic.csv", index=False)
    liv_freq.to_csv(ANALYSIS_DIR / "top1_frequency_liver.csv", index=False)
    comparison.to_csv(ANALYSIS_DIR / "top1_frequency_comparison.csv", index=False)
    write_text(freq_text, ANALYSIS_DIR / "top1_frequency_summary.txt")
    print(freq_text)
    print()

    print("=== Analysis 3: TTK deep-dive ===")
    ttk = analysis_ttk_deepdive()
    write_text(ttk, ANALYSIS_DIR / "ttk_deepdive.txt")
    print(ttk)
    print()

    print("=== Analysis 4: Ranking shifts ===")
    shifts = analysis_ranking_shifts()
    write_text(shifts, ANALYSIS_DIR / "ranking_shifts.txt")
    print(shifts)
    print()

    print(f"All outputs written to: {ANALYSIS_DIR}")


if __name__ == "__main__":
    main()
