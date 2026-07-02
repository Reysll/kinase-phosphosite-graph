"""
Poster Figure 3: LOO prediction performance.
Grouped bar chart comparing node2vec vs spectral median adjusted rank,
across 4 networks. Data are Phase 1/2 multi-kinase results (6,909 trials, 420 candidates).
Lower rank = better prediction (max rank = 420).
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import os

OUT_PATH = "outputs_prediction/analysis/poster_loo_figure.png"

NETWORKS = ["Generic", "Control", "Cancer", "Liver FC"]

NODE2VEC = [188.0, 166.5, 174.0, 165.0]
SPECTRAL = [79.0, 74.0, 74.0, 74.0]

x = np.arange(len(NETWORKS))
width = 0.35

fig, ax = plt.subplots(figsize=(7, 4.5))

COLOR_N2V = "#5B9BD5"  # blue-grey
COLOR_SPE = "#ED7D31"  # orange

bars_n2v = ax.bar(
    x - width / 2,
    NODE2VEC,
    width,
    label="Node2Vec",
    color=COLOR_N2V,
    edgecolor="white",
    linewidth=0.8,
)
bars_spe = ax.bar(
    x + width / 2,
    SPECTRAL,
    width,
    label="Spectral",
    color=COLOR_SPE,
    edgecolor="white",
    linewidth=0.8,
)

for bar, val in zip(bars_n2v, NODE2VEC):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 4,
        f"{val:.0f}",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#333333",
    )

for bar, val in zip(bars_spe, SPECTRAL):
    ax.text(
        bar.get_x() + bar.get_width() / 2,
        bar.get_height() + 4,
        f"{val:.0f}",
        ha="center",
        va="bottom",
        fontsize=8,
        color="#333333",
    )

ax.set_xticks(x)
ax.set_xticklabels(NETWORKS, fontsize=10)
ax.set_ylabel("Median Adjusted Rank", fontsize=10)
ax.set_ylim(0, 230)
ax.set_title(
    "LOO Evaluation: Median Kinase Rank by Network and Embedding\n"
    "(420 candidate kinases; 6,909 multi-kinase trials; lower = better)",
    fontsize=9,
    pad=8,
)
ax.legend(fontsize=9, framealpha=0.9)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.yaxis.set_tick_params(labelsize=9)
ax.axhline(y=210, color="#aaaaaa", linestyle="--", linewidth=0.7, zorder=0)
ax.text(3.62, 212, "random\nbaseline\n(210)", ha="center", va="bottom", fontsize=7, color="#888888")

fig.tight_layout()
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight")
print(f"Wrote: {OUT_PATH}")
