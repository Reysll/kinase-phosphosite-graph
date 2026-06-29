"""
Result 1 (Dr. Ayati, agreed 2026-06-24): average rank per kinase per network,
visualized as a clustergram -- one per embedding strategy (node2vec, spectral,
concat), columns = generic/control/cancer/liver FC.

Consumes the long-format outputs of run_inference_all.py
(outputs_prediction/inference_all/{network}_{embedding}/ranks.csv.gz).

Caveat carried over from run_inference_all.py: the "generic" network only
covers 391/3,228 (12%) of the FC-valid target sites (its graph's site
population is essentially PSP-annotated sites, a small subset of the broader
liver-proteomics-detected set). Generic's average ranks below are computed
over that smaller 391-site subset, not the full 3,228 -- a real asymmetry,
not a bug. See run_inference_all.py module docstring for the full explanation.

IMPROVEMENT (2026-06-25, prompted by "are the heatmaps meaningful?"): the
original absolute-mean-rank clustermap was dominated by the same hub-bias
degeneracy documented in build_rank_ks_test.py -- most kinases' rank barely
varies by site, so their mean rank is close to a fixed per-kinase baseline
regardless of network. Checking outputs_prediction/analysis/rank_heatmap_spectral.csv
confirmed it directly: under spectral, Control/Cancer/Liver FC columns have
near-identical per-kinase values (e.g. TSSK4: 414.0/414.0/414.0) and identical
column-level mean (208.0) and std (~119.9). The clustergram's headline finding
("Cancer+Liver FC cluster tightest, Control next, Generic the outlier") was
largely just restating "Generic is structurally different", not revealing a
liver-cancer-specific signal -- the actual question of interest per Dr.
Ayati's Option 2 framing (email_thread.txt, 5/23: score = score_disease -
score_control).

Two fixes applied:
1. Clustermaps (the *_*.png files) are now restricted to "informative"
   kinases -- same definition as build_rank_ks_test.py's `informative` flag
   (std_rank > MIN_INFORMATIVE_STD in at least one of control/cancer/liver_fc;
   generic excluded from this check because its 391-site sample isn't
   comparable to the other three's full 3,228-site coverage). This removes
   kinases whose mean rank is just a fixed hub baseline from the clustering
   distance calculation. The full CSVs still list ALL kinases (with an
   `informative` column) so nothing is hidden, just not plotted.
2. New delta tables/clustermaps: rank_heatmap_delta_{embedding}.{csv,png}.
   delta = mean_rank[disease network] - mean_rank[Control], for Cancer and
   Liver FC vs Control (Control = healthy raw-abundance correlation network,
   used as baseline per Dr. Ayati's Option 2). Negative delta = kinase ranks
   BETTER (more confidently predicted) under the disease network than under
   Control -- i.e. it "gained predicted relevance". This cancels the
   hub-dominated absolute baseline and isolates the network-dependent shift,
   which is the quantity actually relevant to the liver-cancer question.
   Generic is intentionally excluded from the delta pairs (not just plotting)
   -- its 12%-coverage site population differs from Control's, so a delta
   against it would conflate the coverage asymmetry with a real network
   effect.

Outputs (outputs_prediction/analysis/):
  rank_heatmap_{embedding}.png         -- clustergram, informative kinases x networks, absolute mean rank
  rank_heatmap_{embedding}.csv         -- full average-rank table, ALL kinases, with `informative` column
  rank_heatmap_delta_{embedding}.png   -- clustergram, informative kinases x {Cancer,Liver FC} delta vs Control
  rank_heatmap_delta_{embedding}.csv   -- full delta table, ALL kinases
  rank_heatmap_summary.txt             -- coverage notes + top/bottom kinases (absolute and delta)
"""

from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

from src_prediction.core.config import PRED_OUTPUTS_DIR
from src_prediction.core.io_utils import write_text

INFERENCE_ALL_DIR = PRED_OUTPUTS_DIR / "inference_all"
ANALYSIS_DIR = PRED_OUTPUTS_DIR / "analysis"

NETWORKS = ["generic", "control", "cancer", "liver_fc"]
NETWORK_LABELS = {
    "generic": "Generic",
    "control": "Control",
    "cancer": "Cancer",
    "liver_fc": "Liver FC",
}
EMBEDDINGS = ["node2vec", "spectral", "concat"]

# Same threshold and same full-coverage-only restriction as build_rank_ks_test.py's
# `informative` flag -- see that module's docstring for the hub-bias rationale.
# Generic is excluded here because its 391/3,228 (12%) site coverage isn't
# comparable to the other three networks' full 3,228-site coverage.
MIN_INFORMATIVE_STD = 1.0
INFORMATIVE_NETWORKS = ["control", "cancer", "liver_fc"]

# Per Dr. Ayati's Option 2 framing (email_thread.txt, 5/23): score(kinase) =
# score_disease - score_control. Control = healthy raw-abundance correlation
# network is the baseline; deltas isolate the network-dependent shift instead
# of the hub-dominated absolute rank.
DELTA_BASELINE = "control"
DELTA_TARGETS = ["cancer", "liver_fc"]


def _clean_kinase(node_id: str) -> str:
    return str(node_id).split(":", 1)[-1]


def _load_rank_stats(network: str, embedding: str) -> pd.DataFrame | None:
    """Per-kinase mean and std of rank_in_site for one network/embedding combo."""
    path = INFERENCE_ALL_DIR / f"{network}_{embedding}" / "ranks.csv.gz"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    return df.groupby("kinase_node_id")["rank_in_site"].agg(["mean", "std"])


def build_tables(embedding: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (mean_table, std_table), both indexed by clean kinase name x network."""
    mean_cols, std_cols = {}, {}
    for network in NETWORKS:
        stats = _load_rank_stats(network, embedding)
        if stats is not None:
            mean_cols[NETWORK_LABELS[network]] = stats["mean"]
            std_cols[NETWORK_LABELS[network]] = stats["std"]
    mean_table = pd.DataFrame(mean_cols)
    std_table = pd.DataFrame(std_cols)
    for table in (mean_table, std_table):
        table.index = [_clean_kinase(k) for k in table.index]
        table.index.name = "kinase"
    return mean_table.sort_index(), std_table.sort_index()


def compute_informative_mask(std_table: pd.DataFrame) -> pd.Series:
    cols = [
        NETWORK_LABELS[n] for n in INFORMATIVE_NETWORKS if NETWORK_LABELS[n] in std_table.columns
    ]
    if not cols:
        return pd.Series(True, index=std_table.index)
    return (std_table[cols] > MIN_INFORMATIVE_STD).any(axis=1)


def build_delta_table(mean_table: pd.DataFrame) -> pd.DataFrame:
    baseline_col = NETWORK_LABELS[DELTA_BASELINE]
    if baseline_col not in mean_table.columns:
        return pd.DataFrame(index=mean_table.index)
    deltas = {}
    for target in DELTA_TARGETS:
        target_col = NETWORK_LABELS[target]
        if target_col in mean_table.columns:
            deltas[f"{target_col} - {baseline_col}"] = (
                mean_table[target_col] - mean_table[baseline_col]
            )
    return pd.DataFrame(deltas)


def _save_clustermap(
    table: pd.DataFrame,
    title: str,
    out_path,
    cmap: str = "viridis_r",
    center: float | None = None,
    cbar_label: str = "mean rank_in_site (lower = better)",
    vmin: float | None = None,
    vmax: float | None = None,
) -> tuple[int, int]:
    complete = table.dropna()
    n_total, n_complete = len(table), len(complete)
    if n_complete < 2:
        return n_total, n_complete

    g = sns.clustermap(
        complete,
        cmap=cmap,
        center=center,
        vmin=vmin,
        vmax=vmax,
        figsize=(8, max(10, n_complete * 0.08)),
        cbar_kws={"label": cbar_label},
        yticklabels=n_complete <= 150,
    )
    g.fig.suptitle(title, y=1.02)
    g.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(g.fig)
    print(f"Wrote: {out_path}")
    return n_total, n_complete


def _format_row(kinase: str, row: pd.Series) -> str:
    vals = "  ".join(
        f"{col}={row[col]:.1f}" if pd.notna(row[col]) else f"{col}=NA" for col in row.index
    )
    return f"  {kinase:<14}  {vals}"


def main() -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    summary_lines = [
        "Result 1 -- Average Rank Per Kinase Per Network (clustergram)",
        "=" * 70,
        "",
        "Source: outputs_prediction/inference_all/{network}_{embedding}/ranks.csv.gz",
        "Value: mean rank_in_site across all FC-valid target sites where the",
        "       kinase received a score (own-edge masked for known true KSAs,",
        "       global model otherwise -- see inference_all.py).",
        "",
        "COVERAGE CAVEAT: 'generic' only covers 391/3,228 (12%) of FC-valid",
        "sites -- its graph's sites are essentially the PSP-annotated",
        "population, much smaller than the broader liver-proteomics-detected",
        "set covered by control/cancer/liver FC (100% each). Generic's average",
        "ranks below are over that smaller 391-site subset.",
        "",
        "HUB-BIAS CAVEAT: most kinases' mean rank is close to a fixed",
        "per-kinase baseline regardless of network (same degeneracy as",
        "build_rank_ks_test.py's `informative` flag). Clustermaps below are",
        f"restricted to 'informative' kinases (std_rank > {MIN_INFORMATIVE_STD} in at",
        "least one of control/cancer/liver_fc) -- the full CSVs still list every",
        "kinase, with an `informative` column, so nothing is hidden.",
        "",
        "DELTA TABLES: delta = mean_rank[disease] - mean_rank[Control], for",
        "Cancer and Liver FC vs Control (per Dr. Ayati's Option 2 framing:",
        "score = score_disease - score_control). Negative = kinase ranks",
        "BETTER (more confidently predicted) under the disease network --",
        "i.e. gained predicted relevance. Generic is excluded from delta pairs",
        "(coverage asymmetry would conflate with the network effect).",
        "",
    ]

    for embedding in EMBEDDINGS:
        mean_table, std_table = build_tables(embedding)
        if mean_table.empty:
            summary_lines.append(f"{embedding}: no inference_all outputs found, skipped.\n")
            continue

        informative_mask = compute_informative_mask(std_table)
        n_informative = int(informative_mask.sum())

        out_table = mean_table.copy()
        out_table["informative"] = informative_mask.reindex(out_table.index, fill_value=False)
        csv_path = ANALYSIS_DIR / f"rank_heatmap_{embedding}.csv"
        out_table.to_csv(csv_path)
        print(f"Wrote: {csv_path}")

        plot_table = mean_table[informative_mask.reindex(mean_table.index, fill_value=False)]
        n_total, n_complete = _save_clustermap(
            plot_table,
            f"Average kinase rank per network -- {embedding} embedding (informative kinases only)",
            ANALYSIS_DIR / f"rank_heatmap_{embedding}.png",
        )

        summary_lines += [
            f"--- {embedding} ---",
            f"Kinases with data in >=1 network: {len(mean_table)}",
            f"Informative kinases (std_rank > {MIN_INFORMATIVE_STD}): {n_informative} / {len(mean_table)}",
            f"Informative kinases with complete data in ALL networks (plotted): {n_complete}",
            "",
        ]
        if "Liver FC" in mean_table.columns:
            top10 = (
                mean_table.loc[mean_table.index.isin(plot_table.index)]
                .sort_values("Liver FC")
                .head(10)
            )
            summary_lines.append("Top 10 by lowest mean rank in Liver FC (informative only):")
            for kinase, row in top10.iterrows():
                summary_lines.append(_format_row(kinase, row))
        summary_lines.append("")

        # Delta vs Control
        delta_table = build_delta_table(mean_table)
        if not delta_table.empty:
            out_delta = delta_table.copy()
            out_delta["informative"] = informative_mask.reindex(out_delta.index, fill_value=False)
            delta_csv_path = ANALYSIS_DIR / f"rank_heatmap_delta_{embedding}.csv"
            out_delta.to_csv(delta_csv_path)
            print(f"Wrote: {delta_csv_path}")

            delta_plot_table = delta_table[
                informative_mask.reindex(delta_table.index, fill_value=False)
            ]
            max_abs = delta_plot_table.abs().to_numpy()
            vmax = (
                float(max_abs[~pd.isna(max_abs)].max())
                if max_abs.size and not pd.isna(max_abs).all()
                else None
            )
            _save_clustermap(
                delta_plot_table,
                f"Rank delta vs Control -- {embedding} embedding (informative kinases only)",
                ANALYSIS_DIR / f"rank_heatmap_delta_{embedding}.png",
                cmap="RdBu_r",
                center=0,
                cbar_label="rank delta vs Control (negative = gained relevance in disease network)",
                vmin=-vmax if vmax is not None else None,
                vmax=vmax,
            )

            summary_lines.append(f"Delta vs Control (informative kinases only) -- {embedding}:")
            for col in delta_table.columns:
                ranked = delta_plot_table[col].dropna().sort_values()
                if ranked.empty:
                    continue
                summary_lines.append(f"  Most GAINED relevance ({col}, most negative delta):")
                for kinase, val in ranked.head(10).items():
                    summary_lines.append(f"    {kinase:<14}  delta={val:+.1f}")
                summary_lines.append(f"  Most LOST relevance ({col}, most positive delta):")
                for kinase, val in ranked.tail(10).items():
                    summary_lines.append(f"    {kinase:<14}  delta={val:+.1f}")
            summary_lines.append("")

    write_text("\n".join(summary_lines), ANALYSIS_DIR / "rank_heatmap_summary.txt")
    print(f"Wrote: {ANALYSIS_DIR / 'rank_heatmap_summary.txt'}")


if __name__ == "__main__":
    main()
