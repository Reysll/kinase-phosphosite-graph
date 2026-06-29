"""
Volcano plot for Result 2 KS-test: rank-distribution shift between networks.

X axis: median rank delta (network_A - network_B); negative = ranks better in A
Y axis: -log10(adj_p_value)
Shape: filled = informative (std_rank > 1), hollow = degenerate (flat rank)

Generates 4 plots: both pairs × {spectral, node2vec}.
Spectral is the most biologically interpretable.

Reads: outputs_prediction/analysis/ks_test_{pair}_{embedding}.csv
Writes: outputs_prediction/analysis/ks_volcano_{pair}_{embedding}.png
"""

import pathlib
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ANALYSIS_DIR = pathlib.Path("outputs_prediction/analysis")

NECROPTOSIS = {"IRAK4", "RIPK1", "RIPK3", "TRPM7", "MLKL"}
LABEL_TOP_N = 10  # label top N hits by significance + shift magnitude

PAIRS = [
    ("cancer_vs_control", "Cancer", "Control", "Cancer − Control median rank"),
    ("control_vs_liver_fc", "Control", "Liver FC", "Liver FC − Control median rank"),
]
EMBEDDINGS = ["spectral", "node2vec"]


def _plot_volcano(df: pd.DataFrame, x_label: str, title: str, out_path: pathlib.Path) -> None:
    df = df.copy()
    # X = first_network median - second_network median (use the pair order from CSV cols)
    first_col = [
        c for c in df.columns if c.startswith("median_rank_") and "control" not in c.lower()
    ][0]
    second_col = [c for c in df.columns if c.startswith("median_rank_") and c != first_col][0]
    df["delta"] = df[first_col] - df[second_col]
    df["neg_log_p"] = -np.log10(df["adj_p_value"].clip(lower=1e-300))

    fig, ax = plt.subplots(figsize=(9, 6))

    # degenerate (not informative): small hollow markers
    degen = df[~df["informative"]]
    ax.scatter(
        degen["delta"],
        degen["neg_log_p"],
        color="lightgray",
        alpha=0.3,
        s=12,
        zorder=1,
        label="Non-informative (flat rank)",
    )

    # informative, not significant
    info_ns = df[df["informative"] & (df["adj_p_value"] >= 0.05)]
    ax.scatter(
        info_ns["delta"],
        info_ns["neg_log_p"],
        color="#AAAAAA",
        alpha=0.5,
        s=20,
        zorder=2,
        label="Informative, not significant",
    )

    # informative + significant
    info_sig = df[df["informative"] & (df["adj_p_value"] < 0.05)]
    ax.scatter(
        info_sig["delta"],
        info_sig["neg_log_p"],
        color="#4878CF",
        alpha=0.8,
        s=35,
        zorder=3,
        label="Informative + significant",
    )

    # necroptosis cluster
    necro = df[df["kinase"].isin(NECROPTOSIS)]
    ax.scatter(
        necro["delta"],
        necro["neg_log_p"],
        color="#D62728",
        alpha=0.95,
        s=70,
        zorder=5,
        label="Necroptosis / innate-immune cluster",
    )
    for _, row in necro.iterrows():
        ax.annotate(
            row["kinase"],
            (row["delta"], row["neg_log_p"]),
            textcoords="offset points",
            xytext=(5, 3),
            fontsize=8,
            color="#D62728",
            fontweight="bold",
            clip_on=False,
        )

    # label top hits by |delta| × significance (excluding necroptosis already labeled)
    top_candidates = (
        info_sig[~info_sig["kinase"].isin(NECROPTOSIS)]
        .assign(score=lambda d: d["neg_log_p"] * d["delta"].abs())
        .nlargest(LABEL_TOP_N, "score")
    )
    for _, row in top_candidates.iterrows():
        ax.annotate(
            row["kinase"],
            (row["delta"], row["neg_log_p"]),
            textcoords="offset points",
            xytext=(5, -6),
            fontsize=7.5,
            color="#4878CF",
            clip_on=False,
        )

    # significance threshold line
    sig_line = -np.log10(0.05)
    ax.axhline(sig_line, color="gray", linewidth=0.8, linestyle="--", alpha=0.6)
    ax.text(
        ax.get_xlim()[0], sig_line + 0.05, "adj_p = 0.05", fontsize=7, color="gray", va="bottom"
    )

    ax.axvline(0, color="black", linewidth=0.5, linestyle="--", alpha=0.4)

    ax.set_xlabel(x_label, fontsize=11)
    ax.set_ylabel("−log₁₀(adj p-value)", fontsize=11)
    ax.set_title(title, fontsize=11)
    ax.legend(fontsize=8, loc="upper right", framealpha=0.8)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {out_path}")


def main() -> None:
    for pair_key, net_a, net_b, x_label in PAIRS:
        for emb in EMBEDDINGS:
            csv_path = ANALYSIS_DIR / f"ks_test_{pair_key}_{emb}.csv"
            if not csv_path.exists():
                print(f"Skipping (not found): {csv_path}")
                continue
            df = pd.read_csv(csv_path)
            out_path = ANALYSIS_DIR / f"ks_volcano_{pair_key}_{emb}.png"
            title = (
                f"KS-test: rank distribution shift — {net_a} vs {net_b}\n"
                f"({emb} embedding; informative kinases highlighted)"
            )
            _plot_volcano(df, x_label, title, out_path)


if __name__ == "__main__":
    main()
