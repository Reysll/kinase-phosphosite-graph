"""
Supervisor-facing analyses over all four multi-kinase LOO experiments.

Outputs written to outputs_prediction/analysis/:
  performance_summary.txt         — all four node2vec experiments side-by-side
  top1_frequency_summary.txt      — top-20 most predicted kinases per model
  top1_frequency_comparison.csv   — full side-by-side count table
  per_kinase_topk.txt             — per-kinase top-1/5/10 % across networks
  per_kinase_topk.csv             — machine-readable version
  ttk_deepdive.txt                — TTK-specific accuracy and rank metrics (node2vec + spectral)
  ranking_shifts.txt              — top-20 biggest improvements / worsenings (node2vec)
  phase_comparison.txt            — node2vec vs spectral side-by-side (Phase 1 vs 2)
  category_ranks.txt              — per-category rank breakdown (spectral results)
  ranking_shifts_spectral.txt     — generic vs liver FC under spectral embedding
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src_prediction.core.config import PRED_OUTPUTS_DIR
from src_prediction.core.io_utils import write_text


ANALYSIS_DIR = PRED_OUTPUTS_DIR / "analysis"

# Phase 1 — node2vec experiments
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

# Phase 2 — spectral embedding experiments
SPECTRAL_EXPERIMENT_DIRS = {
    "generic_spectral": PRED_OUTPUTS_DIR / "generic_multi_kinase_spectral",
    "control_spectral": PRED_OUTPUTS_DIR / "control_multi_kinase_spectral",
    "cancer_spectral": PRED_OUTPUTS_DIR / "cancer_multi_kinase_spectral",
    "liver_spectral": PRED_OUTPUTS_DIR / "liver_multi_kinase_spectral",
}

SPECTRAL_EXPERIMENT_LABELS = {
    "generic_spectral": "Generic spectral   (no correlation)",
    "control_spectral": "Control spectral   (ctrl raw abundance corr)",
    "cancer_spectral": "Cancer spectral    (cancer raw abundance corr)",
    "liver_spectral": "Liver FC spectral  (fold-change corr, 80th pct)",
}

SPECTRAL_TAGS = {
    "generic_spectral": "generic_sp",
    "control_spectral": "control_sp",
    "cancer_spectral": "cancer_sp",
    "liver_spectral": "liver_sp",
}

# (node2vec_key, spectral_key, short_label) used in phase comparison table
PHASE_PAIRS = [
    ("generic_multi", "generic_spectral", "Generic"),
    ("control_multi", "control_spectral", "Control"),
    ("cancer_multi", "cancer_spectral", "Cancer"),
    ("liver_multi", "liver_spectral", "Liver FC"),
]

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


def _load_spectral_results(exp_key: str) -> pd.DataFrame | None:
    path = SPECTRAL_EXPERIMENT_DIRS[exp_key] / "results.csv.gz"
    if not path.exists():
        return None
    return pd.read_csv(path)


def _clean_node_id(node_id: str) -> str:
    """PROTEIN:GENE -> GENE,  SITE:GENE-S123 -> GENE-S123"""
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
# Analysis 1: combined performance summary (all four node2vec networks)
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
# Analysis 2: top-1 kinase frequency (all four node2vec networks)
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
# Analysis 3: per-kinase top-k % table across all four node2vec networks
# ---------------------------------------------------------------------------


def analysis_per_kinase_topk(min_trials: int = MIN_TRIALS_PER_KINASE) -> tuple[pd.DataFrame, str]:
    """
    For every kinase with >= min_trials LOO trials across ANY experiment, report
    the percentage of its substrates where it is ranked in top-1, top-5, top-10
    using adjusted_held_out_rank.
    """
    rank_col = "adjusted_held_out_rank"
    ks = [1, 5, 10]

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

    # Sort by liver_fc top-10 descending
    sort_col = "top10_liver_fc" if "top10_liver_fc" in merged.columns else merged.columns[0]
    merged = merged.sort_values(sort_col, ascending=False).reset_index(drop=True)

    tags = [t for t in EXPERIMENT_TAGS.values() if f"n_{t}" in merged.columns]

    lines = [
        f"Per-Kinase Top-k% — adjusted_held_out_rank across all networks (min {min_trials} trials)",
        "=" * 120,
        "Interpretation: what % of a kinase's substrates does the model rank it in top-1 / top-5 / top-10?",
        "  Rank metric: adjusted_held_out_rank (other co-true kinases removed from scored list).",
        "  High % in liver_fc but not generic -> context adds genuine signal for that kinase.",
        "  High % in all networks -> likely degree/hub bias (kinase is well-connected everywhere).",
        "",
    ]

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
# Analysis 4: TTK deep-dive (node2vec and spectral)
# ---------------------------------------------------------------------------


def analysis_ttk_deepdive() -> str:
    lines = [
        "TTK Deep-Dive — LOO Evaluation across all four networks",
        "=" * 70,
        "",
        "TTK (polo-like kinase 4 / Mps1) is a liver-cancer-relevant mitotic kinase.",
        "",
    ]

    def _ttk_section(exp_dict: dict, load_fn) -> list[str]:
        out = []
        for key, label in exp_dict.items():
            df = load_fn(key)
            if df is None:
                out.append(f"{label}  — results not found\n")
                continue

            ttk_rows = df[df["true_kinase_node_id"] == TTK_NODE_ID].copy()
            n = len(ttk_rows)

            out.append(f"{label}")
            out.append("-" * len(label))

            if n == 0:
                out.append("  No TTK trials found in this experiment.\n")
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

            out.append(f"  TTK trials: {n}")
            out.append(_fmt(ho, "held_out_rank     "))
            if adj.notna().any():
                out.append(_fmt(adj, "adjusted_held_out "))
            if bt.notna().any():
                out.append(_fmt(bt, "best_true_rank    "))

            top1_counts = (
                ttk_rows["top1_predicted_kinase"]
                .dropna()
                .apply(_clean_node_id)
                .value_counts()
                .head(10)
            )
            out.append("  Top-10 most predicted as top-1 when TTK is the true kinase:")
            for kinase, count in top1_counts.items():
                marker = " <- TTK correctly predicted" if kinase == "TTK" else ""
                out.append(f"    {kinase:<28}  {count:>4}×{marker}")
            out.append("")
        return out

    lines += ["--- Phase 1: node2vec ---", ""]
    lines += _ttk_section(EXPERIMENT_LABELS, _load_results)
    lines += ["", "--- Phase 2: spectral ---", ""]
    lines += _ttk_section(SPECTRAL_EXPERIMENT_LABELS, _load_spectral_results)

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 5: ranking shift case studies (node2vec generic vs liver FC)
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
        "Ranking Shift Case Studies — Multi-Kinase LOO (node2vec: liver FC vs generic)",
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
# Analysis 6: Phase 1 vs Phase 2 comparison (node2vec vs spectral)
# ---------------------------------------------------------------------------


def analysis_phase_comparison() -> str:
    rank_col = "adjusted_held_out_rank"
    lines = [
        "Phase Comparison — node2vec (Phase 1) vs spectral (Phase 2)",
        "=" * 95,
        "",
        "Metric: adjusted_held_out_rank (other co-true kinases removed before ranking).",
        "",
        f"{'Network':<16}  {'P1 median':>10}  {'P1 mean':>8}  {'P1 top-20':>10}  "
        f"{'P2 median':>10}  {'P2 mean':>8}  {'P2 top-20':>10}  {'d_median':>8}  {'%imp':>6}",
        "-" * 95,
    ]

    for n2v_key, sp_key, label in PHASE_PAIRS:
        df1 = _load_results(n2v_key)
        df2 = _load_spectral_results(sp_key)

        def _stats(df):
            if df is None or rank_col not in df.columns:
                return None, None, None
            col = df[rank_col].dropna()
            n = len(df)
            return col.median(), col.mean(), (col <= 20).sum() / n

        med1, mean1, top20_1 = _stats(df1)
        med2, mean2, top20_2 = _stats(df2)

        if med1 is None or med2 is None:
            lines.append(f"  {label:<14}  — results missing")
            continue

        delta = med1 - med2
        pct_imp = delta / med1 * 100

        lines.append(
            f"  {label:<14}  {med1:>10.1f}  {mean1:>8.1f}  {top20_1:>9.1%}  "
            f"{med2:>10.1f}  {mean2:>8.1f}  {top20_2:>9.1%}  {delta:>+8.1f}  {pct_imp:>5.1f}%"
        )

    lines += [
        "",
        "Notes:",
        "  P1 = Phase 1 (node2vec, DeepWalk mode p=q=1, dimensions=32).",
        "  P2 = Phase 2 (SpectralEmbeddingStrategy, normalized Laplacian, 32 eigenvectors).",
        "  d_median = P1_median - P2_median  (positive = P2 improved).",
        "  %imp = d_median / P1_median x 100.",
    ]
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 7: per-category rank breakdown (spectral results)
# ---------------------------------------------------------------------------


def analysis_category_ranks() -> str:
    cat_rank_col = "adjusted_held_out_rank_in_category"
    global_rank_col = "adjusted_held_out_rank"
    lines = [
        "Per-Category Rank Breakdown — Spectral Embedding (Phase 2)",
        "=" * 80,
        "",
        "Category thresholds: poor <=5 substrates (n~166), average 5-20 (n~133), rich >20 (n~121).",
        "adjusted_held_out_rank_in_category: rank among kinases IN THE SAME CATEGORY only.",
        "adjusted_held_out_rank: global rank (all 420 kinases, co-kinases removed).",
        "",
    ]

    for key, label in SPECTRAL_EXPERIMENT_LABELS.items():
        df = _load_spectral_results(key)
        lines.append(f"{label}")
        lines.append("-" * len(label))

        if df is None:
            lines.append("  — results not found\n")
            continue

        if "held_out_kinase_category" not in df.columns:
            lines.append("  — category columns not present\n")
            continue

        n_total = len(df)
        lines.append(f"  Total trials: {n_total:,}")
        lines.append(
            f"  {'Category':<10}  {'n':>6}  {'global med':>10}  {'cat med':>8}  "
            f"{'cat top-1':>10}  {'cat top-5':>10}  {'within-cat top-1 kinase'}"
        )
        lines.append("  " + "-" * 88)

        for cat in ("poor", "average", "rich"):
            cat_df = df[df["held_out_kinase_category"] == cat]
            n_cat = len(cat_df)
            if n_cat == 0:
                lines.append(f"  {cat:<10}  {0:>6}  — no trials")
                continue

            global_med = (
                cat_df[global_rank_col].median()
                if global_rank_col in cat_df.columns
                else float("nan")
            )
            cat_med = (
                cat_df[cat_rank_col].median() if cat_rank_col in cat_df.columns else float("nan")
            )
            cat_top1_pct = (
                (cat_df[cat_rank_col].dropna() <= 1).sum() / n_cat
                if cat_rank_col in cat_df.columns
                else float("nan")
            )
            cat_top5_pct = (
                (cat_df[cat_rank_col].dropna() <= 5).sum() / n_cat
                if cat_rank_col in cat_df.columns
                else float("nan")
            )

            cat_col = f"top1_{cat}"
            if cat_col in cat_df.columns:
                top1_hub = cat_df[cat_col].dropna().apply(_clean_node_id).value_counts()
                hub_str = (
                    f"{top1_hub.index[0]} {top1_hub.iloc[0] / n_cat:.0%}"
                    if len(top1_hub) > 0
                    else "—"
                )
            else:
                hub_str = "—"

            lines.append(
                f"  {cat:<10}  {n_cat:>6,}  {global_med:>10.1f}  {cat_med:>8.1f}  "
                f"{cat_top1_pct:>9.1%}  {cat_top5_pct:>9.1%}  {hub_str}"
            )
        lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Analysis 8: ranking shifts under spectral (generic vs liver FC)
# ---------------------------------------------------------------------------


def analysis_ranking_shifts_spectral() -> str:
    rank_col = "adjusted_held_out_rank"
    g_df = _load_spectral_results("generic_spectral")
    l_df = _load_spectral_results("liver_spectral")

    if g_df is None or l_df is None:
        return "Spectral results for generic or liver_fc not found."

    g_sub = g_df[
        ["trial_index", "site_node_id", "true_kinase_node_id", rank_col, "top1_predicted_kinase"]
    ].copy()
    l_sub = l_df[["trial_index", rank_col, "top1_predicted_kinase"]].copy()

    g_sub = g_sub.rename(
        columns={rank_col: "generic_adj_rank", "top1_predicted_kinase": "generic_top1"}
    )
    l_sub = l_sub.rename(
        columns={rank_col: "liver_adj_rank", "top1_predicted_kinase": "liver_top1"}
    )

    merged = g_sub.merge(l_sub, on="trial_index", how="inner")
    merged = merged.dropna(subset=["generic_adj_rank", "liver_adj_rank"])
    merged["delta"] = merged["generic_adj_rank"] - merged["liver_adj_rank"]
    merged["liver_improved"] = merged["delta"] > 0

    improved = merged[merged["liver_improved"]].nlargest(20, "delta")
    worsened = merged[~merged["liver_improved"]].nsmallest(20, "delta")

    def fmt_row(row: pd.Series) -> str:
        site = _clean_node_id(str(row["site_node_id"]))
        kinase = _clean_node_id(str(row["true_kinase_node_id"]))
        g_top1 = _clean_node_id(str(row["generic_top1"]))
        l_top1 = _clean_node_id(str(row["liver_top1"]))
        return (
            f"  {site:<28}  true: {kinase:<16}  "
            f"generic: {int(row['generic_adj_rank']):>4}  liver: {int(row['liver_adj_rank']):>4}  "
            f"delta: {row['delta']:>+6.0f}  top-1: {g_top1} -> {l_top1}"
        )

    n_imp = int(merged["liver_improved"].sum())
    n_wors = int((~merged["liver_improved"]).sum())

    lines = [
        "Ranking Shift Case Studies — SPECTRAL: Liver FC vs Generic",
        "=" * 90,
        "",
        "Metric: adjusted_held_out_rank (spectral embedding).",
        "delta = generic_adj_rank - liver_adj_rank  (positive = liver improved).",
        "",
        f"Total trials with both results: {len(merged):,}",
        f"  Improved in liver FC: {n_imp:,}  ({n_imp / len(merged):.1%})",
        f"  Worsened in liver FC: {n_wors:,}  ({n_wors / len(merged):.1%})",
        "",
        "Top-20 BIGGEST IMPROVEMENTS in liver_fc (spectral):",
        "-" * 90,
    ]
    for _, row in improved.iterrows():
        lines.append(fmt_row(row))

    lines += [
        "",
        "Top-20 BIGGEST WORSENINGS in liver_fc (spectral):",
        "-" * 90,
    ]
    for _, row in worsened.iterrows():
        lines.append(fmt_row(row))

    # Per-kinase summary
    kinase_summary = (
        merged.groupby("true_kinase_node_id")
        .agg(
            n_trials=("trial_index", "count"),
            mean_delta=("delta", "mean"),
            n_improved=("liver_improved", "sum"),
        )
        .reset_index()
    )
    kinase_summary["kinase"] = kinase_summary["true_kinase_node_id"].apply(_clean_node_id)
    kinase_summary["pct_improved"] = kinase_summary["n_improved"] / kinase_summary["n_trials"]
    kinase_summary = kinase_summary[kinase_summary["n_trials"] >= MIN_TRIALS_PER_KINASE]

    top_imp = kinase_summary.nlargest(15, "mean_delta")
    top_wors = kinase_summary.nsmallest(15, "mean_delta")

    lines += [
        "",
        f"Kinases MOST IMPROVED in liver_fc spectral (min {MIN_TRIALS_PER_KINASE} trials):",
        f"  {'Kinase':<20}  {'n':>6}  {'mean_delta':>10}  {'n_improved':>10}  {'pct_improved':>12}",
        "  " + "-" * 64,
    ]
    for _, row in top_imp.iterrows():
        lines.append(
            f"  {row['kinase']:<20}  {int(row['n_trials']):>6}  "
            f"{row['mean_delta']:>+10.1f}  {int(row['n_improved']):>10}  {row['pct_improved']:>11.1%}"
        )

    lines += [
        "",
        f"Kinases MOST WORSENED in liver_fc spectral (min {MIN_TRIALS_PER_KINASE} trials):",
        f"  {'Kinase':<20}  {'n':>6}  {'mean_delta':>10}  {'n_improved':>10}  {'pct_improved':>12}",
        "  " + "-" * 64,
    ]
    for _, row in top_wors.iterrows():
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

    print("=== Analysis 1: Performance summary (node2vec, all four networks) ===")
    perf = analysis_performance_summary()
    write_text(perf, ANALYSIS_DIR / "performance_summary.txt")
    print(perf)
    print()

    print("=== Analysis 2: Top-1 kinase frequency (node2vec) ===")
    comparison_df, freq_text = analysis_top1_frequency()
    if not comparison_df.empty:
        comparison_df.to_csv(ANALYSIS_DIR / "top1_frequency_comparison.csv", index=False)
    write_text(freq_text, ANALYSIS_DIR / "top1_frequency_summary.txt")
    print(freq_text)
    print()

    print("=== Analysis 3: Per-kinase top-k% across networks (node2vec) ===")
    topk_df, topk_text = analysis_per_kinase_topk()
    if not topk_df.empty:
        topk_df.to_csv(ANALYSIS_DIR / "per_kinase_topk.csv", index=False)
    write_text(topk_text, ANALYSIS_DIR / "per_kinase_topk.txt")
    print(topk_text)
    print()

    print("=== Analysis 4: TTK deep-dive (node2vec + spectral) ===")
    ttk = analysis_ttk_deepdive()
    write_text(ttk, ANALYSIS_DIR / "ttk_deepdive.txt")
    print(ttk)
    print()

    print("=== Analysis 5: Ranking shifts (node2vec generic vs liver FC) ===")
    shifts = analysis_ranking_shifts()
    write_text(shifts, ANALYSIS_DIR / "ranking_shifts.txt")
    print(shifts)
    print()

    print("=== Analysis 6: Phase comparison (node2vec vs spectral) ===")
    phase_cmp = analysis_phase_comparison()
    write_text(phase_cmp, ANALYSIS_DIR / "phase_comparison.txt")
    print(phase_cmp)
    print()

    print("=== Analysis 7: Per-category rank breakdown (spectral) ===")
    cat_ranks = analysis_category_ranks()
    write_text(cat_ranks, ANALYSIS_DIR / "category_ranks.txt")
    print(cat_ranks)
    print()

    print("=== Analysis 8: Ranking shifts under spectral (generic vs liver FC) ===")
    shifts_sp = analysis_ranking_shifts_spectral()
    write_text(shifts_sp, ANALYSIS_DIR / "ranking_shifts_spectral.txt")
    print(shifts_sp)
    print()

    print(f"All outputs written to: {ANALYSIS_DIR}")


if __name__ == "__main__":
    main()
