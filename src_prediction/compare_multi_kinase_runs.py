from __future__ import annotations

from pathlib import Path
import pandas as pd


OUTPUT_DIR = Path("outputs_prediction")

GENERIC_KEYWORD = "generic_trained_model_multi_kinase"
LIVER_KEYWORD = "liver_trained_model_multi_kinase_fc_corr"


def find_results_file(keyword: str) -> Path:
    matches = sorted(OUTPUT_DIR.rglob(f"*{keyword}*results*.csv.gz"))
    if not matches:
        matches = sorted(OUTPUT_DIR.rglob(f"*{keyword}*.csv.gz"))

    # remove obvious non-results files if possible
    matches = [m for m in matches if "metrics" not in m.name.lower() and "comparison" not in m.name.lower()]

    if not matches:
        raise FileNotFoundError(f"Could not find results file for keyword: {keyword}")

    if len(matches) > 1:
        print(f"Multiple matches found for {keyword}. Using: {matches[0]}")
        for m in matches:
            print(f"  - {m}")

    return matches[0]


def main() -> None:
    generic_path = find_results_file(GENERIC_KEYWORD)
    liver_path = find_results_file(LIVER_KEYWORD)

    print(f"Reading generic results: {generic_path}")
    print(f"Reading liver results:   {liver_path}")

    generic = pd.read_csv(generic_path)
    liver = pd.read_csv(liver_path)

    generic = generic.rename(
        columns={
            "held_out_rank": "generic_held_out_rank",
            "best_true_rank": "generic_best_true_rank",
            "top1_predicted_kinase": "generic_top1",
            "top1_score": "generic_top1_score",
        }
    )

    liver = liver.rename(
        columns={
            "held_out_rank": "liver_held_out_rank",
            "best_true_rank": "liver_best_true_rank",
            "top1_predicted_kinase": "liver_top1",
            "top1_score": "liver_top1_score",
        }
    )

    keep_generic = [
        "fold_index",
        "true_kinase_node_id",
        "site_node_id",
        "held_out_relation",
        "generic_held_out_rank",
        "generic_best_true_rank",
        "n_true_kinases_for_site",
        "all_true_kinases_for_site",
        "generic_top1",
        "generic_top1_score",
        "num_candidate_kinases_scored",
        "site_embedding_present",
        "true_kinase_embedding_present",
    ]

    keep_liver = [
        "fold_index",
        "true_kinase_node_id",
        "site_node_id",
        "held_out_relation",
        "liver_held_out_rank",
        "liver_best_true_rank",
        "n_true_kinases_for_site",
        "all_true_kinases_for_site",
        "liver_top1",
        "liver_top1_score",
        "num_candidate_kinases_scored",
        "site_embedding_present",
        "true_kinase_embedding_present",
    ]

    generic = generic[keep_generic].copy()
    liver = liver[keep_liver].copy()

    compare = generic.merge(
        liver,
        on=["fold_index", "true_kinase_node_id", "site_node_id"],
        how="inner",
        suffixes=("", "_dup"),
    )

    compare["held_out_rank_delta"] = compare["generic_held_out_rank"] - compare["liver_held_out_rank"]
    compare["best_true_rank_delta"] = compare["generic_best_true_rank"] - compare["liver_best_true_rank"]
    compare["top1_changed"] = compare["generic_top1"] != compare["liver_top1"]
    compare["held_out_improved_in_liver"] = compare["held_out_rank_delta"] > 0
    compare["best_true_improved_in_liver"] = compare["best_true_rank_delta"] > 0
    compare["generic_top1_is_known_true"] = compare.apply(
        lambda r: str(r["generic_top1"]) in str(r["all_true_kinases_for_site"]).split(";"),
        axis=1,
    )
    compare["liver_top1_is_known_true"] = compare.apply(
        lambda r: str(r["liver_top1"]) in str(r["all_true_kinases_for_site"]).split(";"),
        axis=1,
    )

    all_out = OUTPUT_DIR / "comparison_generic_vs_liver_multi_kinase_fc_corr.csv.gz"
    changed_out = OUTPUT_DIR / "comparison_changed_top1_only_multi_kinase_fc_corr.csv.gz"
    improved_out = OUTPUT_DIR / "comparison_improved_in_liver_multi_kinase_fc_corr.csv.gz"
    strongest_out = OUTPUT_DIR / "comparison_multi_kinase_known_top_predictions_fc_corr.csv.gz"
    summary_out = OUTPUT_DIR / "comparison_summary_multi_kinase_fc_corr.txt"

    compare.to_csv(all_out, index=False, compression="gzip")
    compare[compare["top1_changed"]].to_csv(changed_out, index=False, compression="gzip")
    compare[compare["held_out_improved_in_liver"]].to_csv(improved_out, index=False, compression="gzip")

    strongest = compare[
        compare["top1_changed"] &
        (compare["generic_top1_is_known_true"] | compare["liver_top1_is_known_true"])
    ].copy()
    strongest = strongest.sort_values(
        ["best_true_rank_delta", "held_out_rank_delta"],
        ascending=[False, False]
    )
    strongest.to_csv(strongest_out, index=False, compression="gzip")

    lines = [
        "Comparison: generic trained model vs liver trained model on multi-kinase sites",
        "=" * 78,
        f"Rows compared: {len(compare)}",
        f"Top1 changed: {int(compare['top1_changed'].sum())}",
        f"Held-out improved in liver: {int(compare['held_out_improved_in_liver'].sum())}",
        f"Best-true improved in liver: {int(compare['best_true_improved_in_liver'].sum())}",
        f"Generic top1 is known true: {int(compare['generic_top1_is_known_true'].sum())}",
        f"Liver top1 is known true: {int(compare['liver_top1_is_known_true'].sum())}",
        f"Mean generic held-out rank: {compare['generic_held_out_rank'].mean():.3f}",
        f"Mean liver held-out rank: {compare['liver_held_out_rank'].mean():.3f}",
        f"Median generic held-out rank: {compare['generic_held_out_rank'].median():.3f}",
        f"Median liver held-out rank: {compare['liver_held_out_rank'].median():.3f}",
        f"Mean generic best-true rank: {compare['generic_best_true_rank'].mean():.3f}",
        f"Mean liver best-true rank: {compare['liver_best_true_rank'].mean():.3f}",
        f"Median generic best-true rank: {compare['generic_best_true_rank'].median():.3f}",
        f"Median liver best-true rank: {compare['liver_best_true_rank'].median():.3f}",
        "",
        f"Wrote: {all_out}",
        f"Wrote: {changed_out}",
        f"Wrote: {improved_out}",
        f"Wrote: {strongest_out}",
        f"Wrote: {summary_out}",
    ]

    summary_text = "\n".join(lines)
    summary_out.write_text(summary_text, encoding="utf-8")

    print()
    print(summary_text)


if __name__ == "__main__":
    main()