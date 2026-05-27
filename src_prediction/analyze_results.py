"""
Supervisor-facing analyses over all four multi-kinase LOO experiments.

Outputs written to outputs_prediction/analysis/:
  performance_summary.txt         — all four experiments side-by-side
  top1_frequency_summary.txt      — top-20 most predicted kinases per model
  top1_frequency_comparison.csv   — full side-by-side count table
  per_kinase_topk.txt             — per-kinase top-1/5/10 % across networks
  per_kinase_topk.csv             — machine-readable version
  ttk_deepdive.txt                — TTK-specific accuracy and rank metrics
  ranking_shifts.txt              — top-20 biggest improvements / worsenings
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src_prediction.config import PRED_OUTPUTS_DIR
from src_prediction.io_utils import write_text


ANALYSIS_DIR = PRED_OUTPUTS_DIR / "analysis"

# Four network experiments compared side-by-side
EXPERIMENT_DIRS = {
    "generic_multi": PRED_OUTPUTS_DIR / "generic_multi_kinase",
    "control_multi": PRED_OUTPUTS_DIR / "control_multi_kinase",
    "cancer_multi": PRED_OUTPUTS_DIR / "cancer_multi_kinase",
    "liver_multi": PRED_OUTPUTS_DIR / "liver_multi_kinase",
}

EXPERIMENT_LABELS = {
    "generic_multi": "Generic            (no correlation)",
    "control_multi": "Control / healthy  (ctrl raw abundance corr)",
    "cancer_multi": "Cancer / tumor     (cancer raw abundance corr)",
    "liver_multi": "Liver FC           (fold-change corr, 80th pct)",
}

# Short column tags for wide tables
EXPERIMENT_TAGS = {
    "generic_multi": "generic",
    "control_multi": "control",
    "cancer_multi": "cancer",
    "liver_multi": "liver_fc",
}

TTK_NODE_ID = "PROTEIN:TTK"
MIN_TRIALS_PER_KINASE = 5  # minimum LOO trials to include a kinase in per-kinase table


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


def _rank_stats(series: pd.Series) -> tuple[str, str]:
    valid = series.dropna()
    if valid.empty:
        return "N/A", "N/A"
    return f"{valid.mean():.1f}", f"{valid.median():.0f}"


# ---------------------------------------------------------------------------
# Analysis 1: combined performance summary (all four networks)
# ---------------------------------------------------------------------------


def analysis_performance_summary() -> str:
    col_w = 50

    hdr = (
        f"{'Experiment':<{col_w}}  {'n':>6}  "
        f"{'Top-1':>7}  {'Top-5':>7}  {'Top-10':>7}  {'Top-20':>7}  "
        f"{'Mean':>8}  {'Median':>8}"
    )
    sep = "-" * len(hdr)

    lines = [
        "Combined Performance Summary — LOO Evaluation (all four networks)",
        "=" * len(hdr),
        "",
    ]

    def _section(rank_col: str, section_title: str) -> list[str]:
        out = ["", sep, section_title, sep]
        for key, label in EXPERIMENT_LABELS.items():
            df = _load_results(key)
            if df is None:
                out.append(f"  {label}  — results not found")
                continue
            if rank_col not in df.columns:
                out.append(f"  {label}  — column '{rank_col}' not present")
                continue
            n = len(df)
            col = df[rank_col]
            mean_r, med_r = _rank_stats(col)
            out.append(
                f"  {label:<{col_w - 2}}"
                f"  {n:>6,}"
                f"  {_top_k_str(col, 1, n):>7}"
                f"  {_top_k_str(col, 5, n):>7}"
                f"  {_top_k_str(col, 10, n):>7}"
                f"  {_top_k_str(col, 20, n):>7}"
                f"  {mean_r:>8}"
                f"  {med_r:>8}"
            )
        return out

    lines.append(hdr)
    lines += _section(
        "held_out_kinase_rank",
        "Held-out rank  (naive: held-out kinase vs all 420 candidates)",
    )
    lines += _section(
        "adjusted_held_out_rank",
        "Adjusted held-out rank  (other co-true kinases removed before ranking — fairer metric)",
    )
    lines += _section(
        "best_true_kinase_rank",
        "Best-true rank  (best-ranked known true kinase for the site)",
    )

    lines += [
        "",
        "Notes:",
        "  candidate pool = 420 kinases",
        "  adjusted_held_out_rank: all OTHER known true kinases for the site are",
        "    removed from the ranked list before computing the held-out kinase rank.",
        "    e.g. if site has k1/k2/k3 and k1 is held out, k2 and k3 are dropped",
        "    from scoring before k1's position is measured.",
        "  best_true_rank: best rank any known true kinase achieves (not adjusted).",
    ]

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 2: top-1 kinase frequency (all four networks)
# ---------------------------------------------------------------------------


def analysis_top1_frequency() -> tuple[pd.DataFrame, str]:
    freq_tables: dict[str, pd.DataFrame] = {}
    for key, tag in EXPERIMENT_TAGS.items():
        df = _load_results(key)
        if df is None:
            continue
        counts = (
            df["top1_predicted_kinase"].dropna().apply(_clean_node_id).value_counts().reset_index()
        )
        counts.columns = ["kinase", f"count_{tag}"]
        counts[f"pct_{tag}"] = counts[f"count_{tag}"] / len(df)
        freq_tables[tag] = counts

    if not freq_tables:
        return pd.DataFrame(), "No results loaded."

    # Merge all into one comparison table
    comparison = None
    for tag, tbl in freq_tables.items():
        if comparison is None:
            comparison = tbl
        else:
            comparison = comparison.merge(tbl, on="kinase", how="outer")
    comparison = comparison.fillna(0)
    for tag in EXPERIMENT_TAGS.values():
        if f"count_{tag}" in comparison.columns:
            comparison[f"count_{tag}"] = comparison[f"count_{tag}"].astype(int)

    # Sort by liver_fc count (primary interest)
    sort_col = "count_liver_fc" if "count_liver_fc" in comparison.columns else comparison.columns[1]
    comparison = comparison.sort_values(sort_col, ascending=False).reset_index(drop=True)

    lines = [
        "Top-1 Predicted Kinase Frequency — Multi-Kinase LOO (per network)",
        "=" * 80,
        "",
    ]

    for key, tag in EXPERIMENT_TAGS.items():
        df = _load_results(key)
        if df is None:
            continue
        n = len(df)
        label = EXPERIMENT_LABELS[key]
        col = f"count_{tag}"
        top20 = freq_tables[tag].head(20)
        lines += [
            f"{label}  (n={n:,})",
            f"  {'Rank':<5}  {'Kinase':<28}  {'Count':>6}  {'%':>6}",
            "  " + "-" * 50,
        ]
        for rank_i, (_, row) in enumerate(top20.iterrows(), 1):
            lines.append(
                f"  {rank_i:<5}  {row['kinase']:<28}  {int(row[col]):>6}  {row[f'pct_{tag}']:>5.1%}"
            )
        lines.append("")

    return comparison, "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 3: per-kinase top-k % table across all four networks
# ---------------------------------------------------------------------------


def analysis_per_kinase_topk(min_trials: int = MIN_TRIALS_PER_KINASE) -> tuple[pd.DataFrame, str]:
    """
    For every kinase with ≥ min_trials LOO trials across ANY experiment, report
    the percentage of its substrates where it is ranked in top-1, top-5, top-10
    using the *adjusted_held_out_rank* (other true co-kinases removed).

    This is the kinase-centric equivalent of the VRK family analysis:
    - High % in both generic and liver -> likely a hub/degree bias, not liver-specific
    - High % in liver but not generic -> genuinely liver-context-boosted
    - High % in generic but not liver -> liver context dilutes the signal
    """
    rank_col = "adjusted_held_out_rank"
    ks = [1, 5, 10]

    all_rows: list[dict] = []

    # Collect per-kinase stats for each network
    per_exp: dict[str, pd.DataFrame] = {}
    for key, tag in EXPERIMENT_TAGS.items():
        df = _load_results(key)
        if df is None or rank_col not in df.columns:
            continue
        grp = df.groupby("true_kinase_node_id")
        stats = grp.apply(
            lambda g: pd.Series(
                {
                    "n_trials": len(g),
                    **{
                        f"top{k}_pct": float((g[rank_col].dropna() <= k).sum()) / len(g) for k in ks
                    },
                    "median_rank": g[rank_col].dropna().median()
                    if g[rank_col].notna().any()
                    else float("nan"),
                }
            )
        ).reset_index()
        stats.columns = ["true_kinase_node_id"] + list(stats.columns[1:])
        stats = stats.rename(
            columns={
                "n_trials": f"n_{tag}",
                **{f"top{k}_pct": f"top{k}_{tag}" for k in ks},
                "median_rank": f"median_{tag}",
            }
        )
        per_exp[tag] = stats

    if not per_exp:
        return pd.DataFrame(), "No results with adjusted_held_out_rank loaded."

    # Merge all experiments
    merged = None
    for tag, stats in per_exp.items():
        if merged is None:
            merged = stats
        else:
            merged = merged.merge(stats, on="true_kinase_node_id", how="outer")

    merged["kinase"] = merged["true_kinase_node_id"].apply(_clean_node_id)
    merged = merged.drop(columns=["true_kinase_node_id"])

    # Keep only kinases with enough trials in at least one experiment
    n_cols = [f"n_{tag}" for tag in EXPERIMENT_TAGS.values() if f"n_{tag}" in merged.columns]
    merged["max_n"] = merged[n_cols].max(axis=1)
    merged = merged[merged["max_n"] >= min_trials].copy()

    # Sort by liver_fc top-10 descending (shows most liver-relevant kinases first)
    sort_col = "top10_liver_fc" if "top10_liver_fc" in merged.columns else merged.columns[0]
    merged = merged.sort_values(sort_col, ascending=False).reset_index(drop=True)

    # Build readable text table
    tags = [t for t in EXPERIMENT_TAGS.values() if f"n_{t}" in merged.columns]

    # Header
    tag_hdrs = "  ".join(f"{'  '.join(['top1', 'top5', 'top10'])} ({t})" for t in tags)
    lines = [
        f"Per-Kinase Top-k% — adjusted_held_out_rank across all networks (min {min_trials} trials)",
        "=" * 120,
        f"Interpretation: what % of a kinase's substrates does the model rank it in top-1 / top-5 / top-10?",
        f"  Rank metric: adjusted_held_out_rank (other co-true kinases removed from scored list).",
        f"  High % in liver_fc but not generic -> context adds genuine signal for that kinase.",
        f"  High % in all networks -> likely degree/hub bias (kinase is well-connected everywhere).",
        "",
    ]

    # Column headers
    hdr_kinase = f"{'Kinase':<22}"
    hdr_n = "  ".join(f"{'n':>4}({t[:4]})" for t in tags)
    hdr_tops = "  ".join(f"{'top1':>5} {'top5':>5} {'top10':>5}  (_{t[:6]}_)" for t in tags)
    lines.append(f"{hdr_kinase}  {hdr_n}  {hdr_tops}")
    lines.append("-" * 160)

    for _, row in merged.iterrows():
        kinase = str(row["kinase"])
        ns = "  ".join(
            f"{int(row[f'n_{t}']):>9}"
            if f"n_{t}" in row.index and pd.notna(row[f"n_{t}"])
            else f"{'—':>9}"
            for t in tags
        )
        tops = "  ".join(
            f"{row[f'top1_{t}']:>5.1%} {row[f'top5_{t}']:>5.1%} {row[f'top10_{t}']:>5.1%}"
            if all(f"top{k}_{t}" in row.index for k in [1, 5, 10])
            else f"{'—':>5} {'—':>5} {'—':>5}"
            for t in tags
        )
        lines.append(f"{kinase:<22}  {ns}  {tops}")

    return merged, "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 4: TTK deep-dive
# ---------------------------------------------------------------------------


def analysis_ttk_deepdive() -> str:
    lines = [
        "TTK Deep-Dive — LOO Evaluation across all four networks",
        "=" * 70,
        "",
        "TTK (polo-like kinase 4 / Mps1) is a liver-cancer-relevant mitotic kinase.",
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
        adj = (
            ttk_rows["adjusted_held_out_rank"]
            if "adjusted_held_out_rank" in ttk_rows.columns
            else pd.Series(dtype=float)
        )
        bt = (
            ttk_rows["best_true_kinase_rank"]
            if "best_true_kinase_rank" in ttk_rows.columns
            else pd.Series(dtype=float)
        )

        def _fmt(series: pd.Series, label_prefix: str) -> str:
            if series.empty or series.isna().all():
                return f"  {label_prefix}: N/A"
            return (
                f"  {label_prefix}: "
                f"top-1 {_top_k_str(series, 1, n)}  "
                f"top-5 {_top_k_str(series, 5, n)}  "
                f"top-10 {_top_k_str(series, 10, n)}  "
                f"mean={series.mean():.1f}  median={series.median():.0f}"
            )

        lines.append(f"  TTK trials: {n}")
        lines.append(_fmt(ho, "held_out_rank     "))
        if adj.notna().any():
            lines.append(_fmt(adj, "adjusted_held_out "))
        if bt.notna().any():
            lines.append(_fmt(bt, "best_true_rank    "))

        top1_counts = (
            ttk_rows["top1_predicted_kinase"].dropna().apply(_clean_node_id).value_counts().head(10)
        )
        lines.append("  Top-10 most predicted as top-1 when TTK is the true kinase:")
        for kinase, count in top1_counts.items():
            marker = " <- TTK correctly predicted" if kinase == "TTK" else ""
            lines.append(f"    {kinase:<28}  {count:>4}×{marker}")
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 5: ranking shift case studies (generic vs liver FC)
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
            if pd.notna(row.get("generic_held_out_rank"))
            else "N/A"
        )
        l_rank = (
            f"{int(row['liver_held_out_rank'])}"
            if pd.notna(row.get("liver_held_out_rank"))
            else "N/A"
        )
        delta = f"{row['held_out_rank_delta']:+.0f}"
        g_top1 = _clean_node_id(str(row.get("generic_top1", "?")))
        l_top1 = _clean_node_id(str(row.get("liver_top1", "?")))
        return (
            f"  {site:<28}  true: {kinase:<16}  "
            f"generic: {g_rank:>4}  liver_fc: {l_rank:>4}  delta: {delta:>6}  "
            f"top-1: {g_top1} -> {l_top1}"
        )

    improved = merged[merged["held_out_rank_delta"] > 0].nlargest(20, "held_out_rank_delta")
    worsened = merged[merged["held_out_rank_delta"] < 0].nsmallest(20, "held_out_rank_delta")

    lines = [
        "Ranking Shift Case Studies — Multi-Kinase LOO (liver FC vs generic)",
        "=" * 90,
        "",
        "delta = generic_rank - liver_rank",
        "  Positive delta -> liver improved ranking (smaller rank = better)",
        "  Negative delta -> liver worsened ranking",
        "",
        f"Top-20 BIGGEST IMPROVEMENTS in liver_fc  (n_improved = {int(merged['held_out_improved_in_liver'].sum()):,}):",
        "-" * 90,
    ]
    for _, row in improved.iterrows():
        lines.append(fmt_row(row))

    lines += [
        "",
        f"Top-20 BIGGEST WORSENINGS in liver_fc  (n_worsened = {int((~merged['held_out_improved_in_liver']).sum()):,}):",
        "-" * 90,
    ]
    for _, row in worsened.iterrows():
        lines.append(fmt_row(row))

    # Per-kinase summary
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
    kinase_summary = kinase_summary[kinase_summary["n_trials"] >= MIN_TRIALS_PER_KINASE]

    top_improved_kinases = kinase_summary.nlargest(15, "mean_delta")
    top_worsened_kinases = kinase_summary.nsmallest(15, "mean_delta")

    lines += [
        "",
        f"Kinases MOST IMPROVED in liver_fc (min {MIN_TRIALS_PER_KINASE} trials, sorted by mean delta):",
        f"  {'Kinase':<20}  {'n':>6}  {'mean_delta':>10}  {'n_improved':>10}  {'pct_improved':>12}",
        "  " + "-" * 64,
    ]
    for _, row in top_improved_kinases.iterrows():
        lines.append(
            f"  {row['kinase']:<20}  {int(row['n_trials']):>6}  "
            f"{row['mean_delta']:>+10.1f}  {int(row['n_improved']):>10}  {row['pct_improved']:>11.1%}"
        )

    lines += [
        "",
        f"Kinases MOST WORSENED in liver_fc (min {MIN_TRIALS_PER_KINASE} trials):",
        f"  {'Kinase':<20}  {'n':>6}  {'mean_delta':>10}  {'n_improved':>10}  {'pct_improved':>12}",
        "  " + "-" * 64,
    ]
    for _, row in top_worsened_kinases.iterrows():
        lines.append(
            f"  {row['kinase']:<20}  {int(row['n_trials']):>6}  "
            f"{row['mean_delta']:>+10.1f}  {int(row['n_improved']):>10}  {row['pct_improved']:>11.1%}"
        )

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Analysis 1: Performance summary (all four networks) ===")
    perf = analysis_performance_summary()
    write_text(perf, ANALYSIS_DIR / "performance_summary.txt")
    print(perf)
    print()

    print("=== Analysis 2: Top-1 kinase frequency ===")
    comparison_df, freq_text = analysis_top1_frequency()
    if not comparison_df.empty:
        comparison_df.to_csv(ANALYSIS_DIR / "top1_frequency_comparison.csv", index=False)
    write_text(freq_text, ANALYSIS_DIR / "top1_frequency_summary.txt")
    print(freq_text)
    print()

    print("=== Analysis 3: Per-kinase top-k% across networks ===")
    topk_df, topk_text = analysis_per_kinase_topk()
    if not topk_df.empty:
        topk_df.to_csv(ANALYSIS_DIR / "per_kinase_topk.csv", index=False)
    write_text(topk_text, ANALYSIS_DIR / "per_kinase_topk.txt")
    print(topk_text)
    print()

    print("=== Analysis 4: TTK deep-dive ===")
    ttk = analysis_ttk_deepdive()
    write_text(ttk, ANALYSIS_DIR / "ttk_deepdive.txt")
    print(ttk)
    print()

    print("=== Analysis 5: Ranking shifts (generic vs liver FC) ===")
    shifts = analysis_ranking_shifts()
    write_text(shifts, ANALYSIS_DIR / "ranking_shifts.txt")
    print(shifts)
    print()

    print(f"All outputs written to: {ANALYSIS_DIR}")


if __name__ == "__main__":
    main()
