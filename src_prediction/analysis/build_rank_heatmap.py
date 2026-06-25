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

Outputs (outputs_prediction/analysis/):
  rank_heatmap_{embedding}.png   -- clustergram, kinases x networks
  rank_heatmap_{embedding}.csv   -- full average-rank table (with NaNs)
  rank_heatmap_summary.txt       -- coverage notes + top/bottom kinases
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


def _clean_kinase(node_id: str) -> str:
    return str(node_id).split(":", 1)[-1]


def _load_avg_rank(network: str, embedding: str) -> pd.Series | None:
    path = INFERENCE_ALL_DIR / f"{network}_{embedding}" / "ranks.csv.gz"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    return df.groupby("kinase_node_id")["rank_in_site"].mean()


def build_table(embedding: str) -> pd.DataFrame:
    cols = {}
    for network in NETWORKS:
        avg = _load_avg_rank(network, embedding)
        if avg is not None:
            cols[NETWORK_LABELS[network]] = avg
    table = pd.DataFrame(cols)
    table.index = [_clean_kinase(k) for k in table.index]
    table.index.name = "kinase"
    return table.sort_index()


def _save_clustermap(table: pd.DataFrame, embedding: str) -> tuple[int, int]:
    complete = table.dropna()
    n_total, n_complete = len(table), len(complete)
    if n_complete < 2:
        return n_total, n_complete

    g = sns.clustermap(
        complete,
        cmap="viridis_r",
        figsize=(8, max(10, n_complete * 0.08)),
        cbar_kws={"label": "mean rank_in_site (lower = better)"},
        yticklabels=n_complete <= 150,
    )
    g.fig.suptitle(f"Average kinase rank per network -- {embedding} embedding", y=1.02)
    out_path = ANALYSIS_DIR / f"rank_heatmap_{embedding}.png"
    g.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(g.fig)
    print(f"Wrote: {out_path}")
    return n_total, n_complete


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
    ]

    for embedding in EMBEDDINGS:
        table = build_table(embedding)
        if table.empty:
            summary_lines.append(f"{embedding}: no inference_all outputs found, skipped.\n")
            continue

        csv_path = ANALYSIS_DIR / f"rank_heatmap_{embedding}.csv"
        table.to_csv(csv_path)
        print(f"Wrote: {csv_path}")

        n_total, n_complete = _save_clustermap(table, embedding)

        summary_lines += [
            f"--- {embedding} ---",
            f"Kinases with data in >=1 network: {n_total}",
            f"Kinases with complete data in ALL networks (used for clustering): {n_complete}",
            "",
        ]
        if "Liver FC" in table.columns:
            top10 = table.sort_values("Liver FC").head(10)
            summary_lines.append("Top 10 by lowest mean rank in Liver FC:")
            for kinase, row in top10.iterrows():
                vals = "  ".join(
                    f"{net}={row[net]:.1f}" if pd.notna(row[net]) else f"{net}=NA"
                    for net in NETWORK_LABELS.values()
                    if net in row.index
                )
                summary_lines.append(f"  {kinase:<14}  {vals}")
        summary_lines.append("")

    write_text("\n".join(summary_lines), ANALYSIS_DIR / "rank_heatmap_summary.txt")
    print(f"Wrote: {ANALYSIS_DIR / 'rank_heatmap_summary.txt'}")


if __name__ == "__main__":
    main()
