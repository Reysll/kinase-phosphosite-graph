"""
Compare generic vs liver model results on the full multi-kinase LOO evaluation
(6,909 trials across 2,459 sites with 2+ known kinases).
"""

from __future__ import annotations

import pandas as pd

from src_prediction.core.config import PRED_OUTPUTS_DIR
from src_prediction.core.io_utils import write_csv_gz, write_text


def run_comparison(generic_path, liver_path, suffix: str, label: str) -> None:
    generic = pd.read_csv(generic_path)
    liver = pd.read_csv(liver_path)

    key_cols = ["trial_index", "true_kinase_node_id", "site_node_id"]

    generic = generic.rename(
        columns={
            "adjusted_held_out_rank": "generic_held_out_rank",
            "best_true_kinase_rank": "generic_best_true_rank",
            "top1_predicted_kinase": "generic_top1",
            "top1_score": "generic_top1_score",
        }
    )

    liver = liver.rename(
        columns={
            "adjusted_held_out_rank": "liver_held_out_rank",
            "best_true_kinase_rank": "liver_best_true_rank",
            "top1_predicted_kinase": "liver_top1",
            "top1_score": "liver_top1_score",
        }
    )

    merged = generic.merge(liver, on=key_cols, how="inner", suffixes=("", "_dup"))

    merged["held_out_rank_delta"] = merged["generic_held_out_rank"] - merged["liver_held_out_rank"]
    merged["best_true_rank_delta"] = (
        merged["generic_best_true_rank"] - merged["liver_best_true_rank"]
    )
    merged["top1_changed"] = merged["generic_top1"] != merged["liver_top1"]
    merged["held_out_improved_in_liver"] = (
        merged["liver_held_out_rank"] < merged["generic_held_out_rank"]
    )
    merged["best_true_improved_in_liver"] = (
        merged["liver_best_true_rank"] < merged["generic_best_true_rank"]
    )

    comparison_out = PRED_OUTPUTS_DIR / f"comparison_multi_kinase_generic_vs_liver{suffix}.csv.gz"
    changed_top1_out = PRED_OUTPUTS_DIR / f"comparison_multi_kinase_changed_top1{suffix}.csv.gz"
    improved_out = PRED_OUTPUTS_DIR / f"comparison_multi_kinase_improved_in_liver{suffix}.csv.gz"
    summary_out = PRED_OUTPUTS_DIR / f"comparison_multi_kinase_summary{suffix}.txt"

    changed_top1 = merged.loc[merged["top1_changed"]].copy()
    improved = (
        merged.loc[merged["held_out_improved_in_liver"]]
        .sort_values("held_out_rank_delta", ascending=False)
        .copy()
    )
    worsened = (
        merged.loc[~merged["held_out_improved_in_liver"]]
        .sort_values("held_out_rank_delta", ascending=True)
        .copy()
    )

    n_trials = len(merged)
    n_candidates = (
        merged["num_candidate_kinases_scored"].max()
        if "num_candidate_kinases_scored" in merged.columns
        else "?"
    )

    def top_k(col: str, k: int) -> str:
        valid = merged[col].dropna()
        n = len(valid)
        return f"{(valid <= k).mean():.1%} ({int((valid <= k).sum())}/{n})"

    summary_lines = [
        f"Comparison: generic model vs liver model — full multi-kinase LOO ({label}, 80% fold-change site correlation)",
        "=" * 100,
        f"LOO trials compared:             {n_trials:,}",
        f"Candidate kinases:               {n_candidates}",
        "",
        "Held-out kinase rank (the removed kinase-site edge)",
        f"  Top-1  generic / liver:        {top_k('generic_held_out_rank', 1)}  /  {top_k('liver_held_out_rank', 1)}",
        f"  Top-5  generic / liver:        {top_k('generic_held_out_rank', 5)}  /  {top_k('liver_held_out_rank', 5)}",
        f"  Top-10 generic / liver:        {top_k('generic_held_out_rank', 10)}  /  {top_k('liver_held_out_rank', 10)}",
        f"  Top-20 generic / liver:        {top_k('generic_held_out_rank', 20)}  /  {top_k('liver_held_out_rank', 20)}",
        f"  Mean rank  generic / liver:    {merged['generic_held_out_rank'].mean():.1f}  /  {merged['liver_held_out_rank'].mean():.1f}",
        f"  Median rank generic / liver:   {merged['generic_held_out_rank'].median():.1f}  /  {merged['liver_held_out_rank'].median():.1f}",
        "",
        "Best-true kinase rank (best-ranked among all known kinases for the site)",
        f"  Top-1  generic / liver:        {top_k('generic_best_true_rank', 1)}  /  {top_k('liver_best_true_rank', 1)}",
        f"  Top-5  generic / liver:        {top_k('generic_best_true_rank', 5)}  /  {top_k('liver_best_true_rank', 5)}",
        f"  Top-10 generic / liver:        {top_k('generic_best_true_rank', 10)}  /  {top_k('liver_best_true_rank', 10)}",
        f"  Top-20 generic / liver:        {top_k('generic_best_true_rank', 20)}  /  {top_k('liver_best_true_rank', 20)}",
        f"  Mean rank  generic / liver:    {merged['generic_best_true_rank'].mean():.1f}  /  {merged['liver_best_true_rank'].mean():.1f}",
        f"  Median rank generic / liver:   {merged['generic_best_true_rank'].median():.1f}  /  {merged['liver_best_true_rank'].median():.1f}",
        "",
        "Per-trial changes",
        f"  Top-1 prediction changed:          {int(merged['top1_changed'].sum()):,} / {n_trials:,}  ({merged['top1_changed'].mean():.1%})",
        f"  Held-out rank improved (liver):    {int(merged['held_out_improved_in_liver'].sum()):,} / {n_trials:,}  ({merged['held_out_improved_in_liver'].mean():.1%})",
        f"  Best-true rank improved (liver):   {int(merged['best_true_improved_in_liver'].sum()):,} / {n_trials:,}  ({merged['best_true_improved_in_liver'].mean():.1%})",
        "",
        "Mean held-out rank delta (generic - liver, positive = liver improved):",
        f"  All trials:                    {merged['held_out_rank_delta'].mean():.2f}",
        f"  Trials improved in liver:      {improved['held_out_rank_delta'].mean():.2f}  (n={len(improved):,})",
        f"  Trials worsened in liver:      {worsened['held_out_rank_delta'].mean():.2f}  (n={len(worsened):,})",
    ]
    summary_text = "\n".join(summary_lines)

    write_csv_gz(merged, comparison_out)
    write_csv_gz(changed_top1, changed_top1_out)
    write_csv_gz(improved, improved_out)
    write_text(summary_text, summary_out)

    print(summary_text)
    print()
    print(f"Wrote: {comparison_out}")
    print(f"Wrote: {changed_top1_out}")
    print(f"Wrote: {improved_out}")
    print(f"Wrote: {summary_out}")


def main() -> None:
    run_comparison(
        generic_path=PRED_OUTPUTS_DIR / "generic_multi_kinase" / "results.csv.gz",
        liver_path=PRED_OUTPUTS_DIR / "liver_multi_kinase" / "results.csv.gz",
        suffix="",
        label="node2vec",
    )
    print()
    run_comparison(
        generic_path=PRED_OUTPUTS_DIR / "generic_multi_kinase_spectral" / "results.csv.gz",
        liver_path=PRED_OUTPUTS_DIR / "liver_multi_kinase_spectral" / "results.csv.gz",
        suffix="_spectral",
        label="spectral",
    )


if __name__ == "__main__":
    main()
