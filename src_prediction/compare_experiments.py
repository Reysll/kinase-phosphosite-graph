from __future__ import annotations

import pandas as pd

from src_prediction.config import PRED_OUTPUTS_DIR
from src_prediction.io_utils import write_csv_gz, write_text


def main() -> None:
    generic_results_path = PRED_OUTPUTS_DIR / "generic_trained_model_debug" / "results.csv.gz"
    liver_results_path = PRED_OUTPUTS_DIR / "liver_trained_model_fc_corr_debug" / "results.csv.gz"

    generic = pd.read_csv(generic_results_path)
    liver = pd.read_csv(liver_results_path)

    key_cols = ["trial_index", "true_kinase_node_id", "site_node_id"]

    generic = generic.rename(
        columns={
            "held_out_kinase_rank": "generic_held_out_rank",
            "best_true_kinase_rank": "generic_best_true_rank",
            "top1_predicted_kinase": "generic_top1",
            "top1_score": "generic_top1_score",
        }
    )

    liver = liver.rename(
        columns={
            "held_out_kinase_rank": "liver_held_out_rank",
            "best_true_kinase_rank": "liver_best_true_rank",
            "top1_predicted_kinase": "liver_top1",
            "top1_score": "liver_top1_score",
        }
    )

    merged = generic.merge(
        liver,
        on=key_cols,
        how="inner",
        suffixes=("", "_dup"),
    )

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

    comparison_out = PRED_OUTPUTS_DIR / "comparison_generic_vs_liver_fc_corr.csv.gz"
    changed_top1_out = PRED_OUTPUTS_DIR / "comparison_changed_top1_only_fc_corr.csv.gz"
    improved_out = PRED_OUTPUTS_DIR / "comparison_improved_in_liver_fc_corr.csv.gz"
    summary_out = PRED_OUTPUTS_DIR / "comparison_summary_fc_corr.txt"

    changed_top1 = merged.loc[merged["top1_changed"]].copy()
    improved = (
        merged.loc[merged["held_out_improved_in_liver"]]
        .sort_values("held_out_rank_delta", ascending=False)
        .copy()
    )

    n_trials = len(merged)
    n_candidates = (
        merged["num_candidate_kinases_scored_x"].max()
        if "num_candidate_kinases_scored_x" in merged.columns
        else "?"
    )

    def top_k(col, k):
        return f"{(merged[col] <= k).mean():.1%} ({int((merged[col] <= k).sum())}/{n_trials})"

    summary_lines = [
        "Comparison: generic trained model vs liver trained model (80% fold-change site correlation)",
        "============================================================================================",
        f"LOO trials compared:             {n_trials}",
        f"Candidate kinases:               {n_candidates}",
        "",
        "Held-out kinase rank",
        f"  Top-1  generic / liver:        {top_k('generic_held_out_rank', 1)}  /  {top_k('liver_held_out_rank', 1)}",
        f"  Top-5  generic / liver:        {top_k('generic_held_out_rank', 5)}  /  {top_k('liver_held_out_rank', 5)}",
        f"  Top-10 generic / liver:        {top_k('generic_held_out_rank', 10)}  /  {top_k('liver_held_out_rank', 10)}",
        f"  Top-20 generic / liver:        {top_k('generic_held_out_rank', 20)}  /  {top_k('liver_held_out_rank', 20)}",
        f"  Mean rank  generic / liver:    {merged['generic_held_out_rank'].mean():.1f}  /  {merged['liver_held_out_rank'].mean():.1f}",
        f"  Median rank generic / liver:   {merged['generic_held_out_rank'].median():.1f}  /  {merged['liver_held_out_rank'].median():.1f}",
        "",
        "Best true kinase rank (multi-kinase sites)",
        f"  Top-1  generic / liver:        {top_k('generic_best_true_rank', 1)}  /  {top_k('liver_best_true_rank', 1)}",
        f"  Top-5  generic / liver:        {top_k('generic_best_true_rank', 5)}  /  {top_k('liver_best_true_rank', 5)}",
        f"  Top-10 generic / liver:        {top_k('generic_best_true_rank', 10)}  /  {top_k('liver_best_true_rank', 10)}",
        f"  Mean rank  generic / liver:    {merged['generic_best_true_rank'].mean():.1f}  /  {merged['liver_best_true_rank'].mean():.1f}",
        f"  Median rank generic / liver:   {merged['generic_best_true_rank'].median():.1f}  /  {merged['liver_best_true_rank'].median():.1f}",
        "",
        "Per-trial changes",
        f"  Top-1 prediction changed:      {int(merged['top1_changed'].sum())} / {n_trials}",
        f"  Held-out rank improved (liver): {int(merged['held_out_improved_in_liver'].sum())} / {n_trials}",
        f"  Best-true rank improved (liver):{int(merged['best_true_improved_in_liver'].sum())} / {n_trials}",
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


if __name__ == "__main__":
    main()
