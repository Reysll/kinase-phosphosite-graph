from __future__ import annotations

import pandas as pd

from src_prediction.config import PRED_OUTPUTS_DIR
from src_prediction.io_utils import write_csv_gz, write_text


def main() -> None:
    generic_results_path = PRED_OUTPUTS_DIR / "generic_cosine_parallel_baseline" / "results.csv.gz"
    liver_results_path = PRED_OUTPUTS_DIR / "liver_cosine_with_site_corr" / "results.csv.gz"

    generic = pd.read_csv(generic_results_path)
    liver = pd.read_csv(liver_results_path)

    key_cols = ["fold_index", "true_kinase_node_id", "site_node_id"]

    generic = generic.rename(
        columns={
            "rank_of_true_kinase": "generic_rank",
            "top1_predicted_kinase": "generic_top1",
            "top1_score": "generic_top1_score",
        }
    )

    liver = liver.rename(
        columns={
            "rank_of_true_kinase": "liver_rank",
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

    merged["rank_delta"] = merged["generic_rank"] - merged["liver_rank"]
    merged["top1_changed"] = merged["generic_top1"] != merged["liver_top1"]
    merged["improved_in_liver"] = merged["liver_rank"] < merged["generic_rank"]
    merged["worsened_in_liver"] = merged["liver_rank"] > merged["generic_rank"]

    comparison_out = PRED_OUTPUTS_DIR / "comparison_generic_vs_liver_site_corr.csv.gz"
    changed_top1_out = PRED_OUTPUTS_DIR / "comparison_changed_top1_only.csv.gz"
    improved_out = PRED_OUTPUTS_DIR / "comparison_improved_in_liver.csv.gz"
    summary_out = PRED_OUTPUTS_DIR / "comparison_summary.txt"

    changed_top1 = merged.loc[merged["top1_changed"]].copy()
    improved = merged.loc[merged["improved_in_liver"]].sort_values("rank_delta", ascending=False).copy()

    summary_lines = [
        "Comparison: generic vs liver + site correlation",
        "=============================================",
        f"Rows compared: {len(merged):,}",
        f"Top1 changed: {int(merged['top1_changed'].sum()):,}",
        f"Improved in liver: {int(merged['improved_in_liver'].sum()):,}",
        f"Worsened in liver: {int(merged['worsened_in_liver'].sum()):,}",
        f"Mean generic rank: {merged['generic_rank'].mean():.3f}",
        f"Mean liver rank: {merged['liver_rank'].mean():.3f}",
        f"Median generic rank: {merged['generic_rank'].median():.3f}",
        f"Median liver rank: {merged['liver_rank'].median():.3f}",
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