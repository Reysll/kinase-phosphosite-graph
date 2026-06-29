"""
Spotlight scatter plot: predicted-relevance shift for disease vs. healthy baseline.

X axis: delta rank (Cancer − Control, spectral)
Y axis: delta rank (Liver FC − Control, spectral)
Negative on both axes = kinase gained predicted relevance in BOTH disease contexts.

Reads: outputs_prediction/analysis/rank_heatmap_delta_spectral.csv
Writes: outputs_prediction/analysis/spotlight_necroptosis.png
"""

import pathlib
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

DELTA_CSV = pathlib.Path("outputs_prediction/analysis/rank_heatmap_delta_spectral.csv")
OUT_PNG = pathlib.Path("outputs_prediction/analysis/spotlight_necroptosis.png")

# Kinase groups for the story
NECROPTOSIS = {"IRAK4", "RIPK1", "RIPK3", "TRPM7", "MLKL"}
# CDK1/PIK3R1 dominate all sites equally → near-zero delta, informative=False
# Draw from full df so they appear at the origin as a reference anchor
HUBS = {"CDK1", "PIK3R1"}
# KSEA-significant kinases that barely move in LOO (biologically active but structurally buried)
KSEA_NOTABLE = {"PRKACA", "MAPK12"}

ALL_LABELED = NECROPTOSIS | HUBS | KSEA_NOTABLE


def main() -> None:
    df = pd.read_csv(DELTA_CSV)
    informative = df[df["informative"]].copy()

    x_col = "Cancer - Control"
    y_col = "Liver FC - Control"

    fig, ax = plt.subplots(figsize=(8, 7))

    # --- background: all informative kinases not in any labeled group ---
    bg = informative[~informative["kinase"].isin(ALL_LABELED)]
    ax.scatter(
        bg[x_col],
        bg[y_col],
        color="lightgray",
        alpha=0.35,
        s=18,
        zorder=1,
    )

    def _scatter_group(subset, color, marker, size, zorder, offset=(5, 3), font_size=8.5):
        if subset.empty:
            return
        ax.scatter(
            subset[x_col],
            subset[y_col],
            color=color,
            alpha=0.9,
            s=size,
            zorder=zorder,
            marker=marker,
        )
        for _, row in subset.iterrows():
            ax.annotate(
                row["kinase"],
                (row[x_col], row[y_col]),
                textcoords="offset points",
                xytext=offset,
                fontsize=font_size,
                color=color,
                fontweight="bold",
                clip_on=False,
            )

    # Hub kinases drawn from FULL df (they are informative=False, near origin)
    _scatter_group(
        df[df["kinase"].isin(HUBS)], color="#4878CF", marker="D", size=65, zorder=3, offset=(5, 4)
    )

    # KSEA-notable kinases from informative df
    _scatter_group(
        informative[informative["kinase"].isin(KSEA_NOTABLE)],
        color="#D4A800",
        marker="s",
        size=60,
        zorder=3,
        font_size=8,
    )

    # Necroptosis cluster — the main finding
    _scatter_group(
        informative[informative["kinase"].isin(NECROPTOSIS)],
        color="#D62728",
        marker="o",
        size=100,
        zorder=4,
    )

    # reference lines at zero
    ax.axhline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)
    ax.axvline(0, color="black", linewidth=0.6, linestyle="--", alpha=0.5)

    # bottom-left quadrant annotation
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_ylim()
    ax.text(
        xmin + 0.02 * (xmax - xmin),
        ymin + 0.02 * (ymax - ymin),
        "Gained relevance in\nboth disease contexts",
        fontsize=7.5,
        color="#D62728",
        va="bottom",
        ha="left",
        style="italic",
    )

    ax.set_xlabel("ΔRank: Cancer − Control (spectral embedding)", fontsize=11)
    ax.set_ylabel("ΔRank: Liver FC − Control (spectral embedding)", fontsize=11)
    ax.set_title(
        "Kinase predicted-relevance shift: disease vs. healthy baseline\n"
        "(negative = gained relevance; spectral embedding; informative kinases + hub anchors)",
        fontsize=10,
    )

    legend_handles = [
        mpatches.Patch(
            color="#D62728", label="Necroptosis / innate-immune (IRAK4, RIPK1, RIPK3, TRPM7, MLKL)"
        ),
        mpatches.Patch(
            color="#4878CF",
            label="Hub kinases CDK1, PIK3R1 — degenerate rank, expected near origin",
        ),
        mpatches.Patch(
            color="#D4A800",
            label="KSEA-significant kinases (biologically active, barely shift in LOO)",
        ),
        mpatches.Patch(
            color="lightgray", label=f"All other informative kinases (n={len(informative)})"
        ),
    ]
    ax.legend(
        handles=legend_handles, fontsize=7.5, loc="upper right", framealpha=0.85, borderpad=0.7
    )

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {OUT_PNG}")

    # Print values for the highlighted kinases
    highlighted = df[df["kinase"].isin(ALL_LABELED)].copy()
    highlighted = highlighted.sort_values("Cancer - Control")
    print("\nDelta values for highlighted kinases:")
    print(
        highlighted[["kinase", "Cancer - Control", "Liver FC - Control", "informative"]].to_string(
            index=False
        )
    )


if __name__ == "__main__":
    main()
