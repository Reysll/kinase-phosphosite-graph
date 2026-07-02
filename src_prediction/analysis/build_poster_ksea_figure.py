"""
Poster Figure 5: KSEA enrichment scores for liver FC dataset.
Horizontal bar chart of all kinases with adj_p < 0.05 (Benjamini-Hochberg).
Positive scores = substrates are hyper-phosphorylated in tumor vs. control.
Negative scores = substrates are hypo-phosphorylated.
"""

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import os

KSEA_CSV = "outputs_prediction/ksea/ksea_psp.csv"
OUT_PATH = "outputs_prediction/analysis/poster_ksea_figure.png"

df = pd.read_csv(KSEA_CSV)
sig = df[df["adj_p_value"] < 0.05].copy()
sig = sig.sort_values("ksea_score", ascending=True)

colors = ["#ED7D31" if s > 0 else "#5B9BD5" for s in sig["ksea_score"]]

fig, ax = plt.subplots(figsize=(6, max(3.5, len(sig) * 0.32)))

bars = ax.barh(
    sig["kinase"], sig["ksea_score"], color=colors, edgecolor="white", linewidth=0.6, height=0.7
)

for bar, val, adj_p in zip(bars, sig["ksea_score"], sig["adj_p_value"]):
    if adj_p < 0.001:
        stars = "***"
    elif adj_p < 0.01:
        stars = "**"
    else:
        stars = "*"
    x_pos = val + (0.12 if val >= 0 else -0.12)
    ha = "left" if val >= 0 else "right"
    ax.text(
        x_pos,
        bar.get_y() + bar.get_height() / 2,
        stars,
        ha=ha,
        va="center",
        fontsize=7,
        color="#333333",
    )

ax.axvline(x=0, color="#444444", linewidth=0.8)
ax.set_xlabel("KSEA Score", fontsize=10)
ax.set_title(
    "Kinase Activity Enrichment (KSEA) in Liver Cancer\n"
    "Significant kinases only (adj. p < 0.05, BH correction; n=17/64 scoreable)",
    fontsize=9,
    pad=8,
)
ax.tick_params(axis="y", labelsize=8)
ax.tick_params(axis="x", labelsize=8)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

from matplotlib.patches import Patch

legend_elements = [
    Patch(facecolor="#ED7D31", label="Up-regulated in tumor"),
    Patch(facecolor="#5B9BD5", label="Down-regulated in tumor"),
]
ax.legend(handles=legend_elements, fontsize=8, framealpha=0.9, loc="lower right")

ax.text(
    0.99,
    -0.13,
    "* p<0.05   ** p<0.01   *** p<0.001",
    transform=ax.transAxes,
    ha="right",
    va="top",
    fontsize=7,
    color="#666666",
)

fig.tight_layout()
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
fig.savefig(OUT_PATH, dpi=300, bbox_inches="tight")
print(f"Wrote: {OUT_PATH}")
