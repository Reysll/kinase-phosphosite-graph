"""
KSEA vs LOO rank scatter: illustrates the hub-bias disconnect.

X axis: KSEA score (Wiredja 2017) from PSP known substrates
Y axis: adjusted median LOO rank under spectral liver-FC network

Kinases with strong KSEA signal but high LOO ranks are invisible to the model
due to hub bias -- they are biologically active but structurally buried.

Reads:
  outputs_prediction/ksea/ksea_psp.csv
  outputs_prediction/liver_multi_kinase_spectral/results.csv.gz
Writes:
  outputs_prediction/analysis/ksea_vs_loo_rank.png
"""

import pathlib
import pandas as pd
import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

KSEA_CSV = pathlib.Path("outputs_prediction/ksea/ksea_psp.csv")
LOO_CSV = pathlib.Path("outputs_prediction/liver_multi_kinase_spectral/results.csv.gz")
OUT_PNG = pathlib.Path("outputs_prediction/analysis/ksea_vs_loo_rank.png")

NECROPTOSIS = {"IRAK4", "RIPK1", "RIPK3", "TRPM7", "MLKL"}
MIN_LOO_TRIALS = 5  # require at least this many LOO trials to report a median rank


def main() -> None:
    ksea = pd.read_csv(KSEA_CSV)
    loo = pd.read_csv(LOO_CSV)

    # per-kinase adjusted median rank from LOO results
    loo["kinase"] = loo["true_kinase_node_id"].str.replace("PROTEIN:", "", regex=False)
    loo_rank = (
        loo.groupby("kinase")
        .agg(
            median_adj_rank=("adjusted_held_out_rank", "median"),
            n_trials=("adjusted_held_out_rank", "count"),
        )
        .reset_index()
    )
    loo_rank = loo_rank[loo_rank["n_trials"] >= MIN_LOO_TRIALS]

    # join on kinase name
    merged = ksea.merge(loo_rank, on="kinase", how="inner")
    if merged.empty:
        print("No rows after merge — check kinase name alignment between KSEA and LOO CSVs.")
        return

    sig = merged["adj_p_value"] < 0.05
    necro = merged["kinase"].isin(NECROPTOSIS)

    fig, ax = plt.subplots(figsize=(9, 6))

    # not significant
    ax.scatter(
        merged.loc[~sig, "ksea_score"],
        merged.loc[~sig, "median_adj_rank"],
        color="lightgray",
        alpha=0.5,
        s=20,
        zorder=1,
        label="KSEA not significant",
    )

    # significant
    ax.scatter(
        merged.loc[sig, "ksea_score"],
        merged.loc[sig, "median_adj_rank"],
        color="#D62728",
        alpha=0.85,
        s=55,
        zorder=3,
        label="KSEA adj_p < 0.05",
    )

    # necroptosis (may not be in KSEA if insufficient substrates)
    if necro.any():
        ax.scatter(
            merged.loc[necro, "ksea_score"],
            merged.loc[necro, "median_adj_rank"],
            color="darkorange",
            alpha=0.95,
            s=80,
            marker="^",
            zorder=4,
            label="Necroptosis cluster",
        )
        for _, row in merged[necro].iterrows():
            ax.annotate(
                row["kinase"],
                (row["ksea_score"], row["median_adj_rank"]),
                textcoords="offset points",
                xytext=(5, 3),
                fontsize=8,
                color="darkorange",
                fontweight="bold",
                clip_on=False,
            )

    # label significant kinases
    for _, row in merged[sig & ~necro].iterrows():
        ax.annotate(
            row["kinase"],
            (row["ksea_score"], row["median_adj_rank"]),
            textcoords="offset points",
            xytext=(5, -5),
            fontsize=7.5,
            color="#D62728",
            clip_on=False,
        )

    # reference line at random baseline (~210/420)
    ax.axhline(210, color="gray", linewidth=0.8, linestyle="--", alpha=0.5)
    ax.text(
        ax.get_xlim()[0], 212, "random baseline (~210)", fontsize=7.5, color="gray", va="bottom"
    )

    ax.set_xlabel("KSEA score (Wiredja 2017; PSP substrates in liver dataset)", fontsize=11)
    ax.set_ylabel("Adjusted median LOO rank — spectral liver FC (lower = better)", fontsize=11)
    ax.set_title(
        "KSEA-active kinases vs. predicted LOO rank\n"
        "(hub bias: biologically active kinases are not necessarily top-ranked)",
        fontsize=11,
    )
    ax.legend(fontsize=9, loc="upper right", framealpha=0.8)
    ax.invert_yaxis()  # rank 1 at top

    plt.tight_layout()
    plt.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"Saved {OUT_PNG}")
    print(f"\nMerged {len(merged)} kinases (KSEA + LOO overlap).")
    print(f"KSEA significant (adj_p<0.05): {sig.sum()} kinases")
    print(
        merged[sig][["kinase", "ksea_score", "median_adj_rank", "adj_p_value"]]
        .sort_values("ksea_score", ascending=False)
        .to_string(index=False)
    )


if __name__ == "__main__":
    main()
