# src_prediction/check_ttk_focus.py

from pathlib import Path
import pandas as pd


OUTPUT_DIR = Path("outputs_prediction")
GENERIC_RESULTS = OUTPUT_DIR / "generic_trained_model_multi_kinase_results.csv.gz"
LIVER_RESULTS = OUTPUT_DIR / "liver_trained_model_multi_kinase_fc_corr_results.csv.gz"

# fallback comparison file if you already generated it
COMPARE_FILE = OUTPUT_DIR / "comparison_generic_vs_liver_multi_kinase_fc_corr.csv.gz"

TTK = "PROTEIN:TTK"


def load_first_existing(paths):
    for p in paths:
        if p.exists():
            return pd.read_csv(p)
    raise FileNotFoundError(f"None of these files were found: {paths}")


def main():
    generic = load_first_existing([
        GENERIC_RESULTS,
        OUTPUT_DIR / "generic_trained_model_multi_kinase_results.csv",
    ])

    liver = load_first_existing([
        LIVER_RESULTS,
        OUTPUT_DIR / "liver_trained_model_multi_kinase_fc_corr_results.csv",
    ])

    # 1. How many known TTK substrates are in the evaluated graph/result set
    ttk_true = liver[liver["true_kinase_node_id"] == TTK].copy()
    n_ttk_true_rows = len(ttk_true)
    n_ttk_true_sites = ttk_true["site_node_id"].nunique()

    # 2. In how many liver predictions is TTK top1
    liver_ttk_top1 = liver[liver["top1_predicted_kinase"] == TTK].copy()
    n_liver_ttk_top1_rows = len(liver_ttk_top1)
    n_liver_ttk_top1_sites = liver_ttk_top1["site_node_id"].nunique()

    # 3. Is it also top1 in generic for those same rows/sites
    merged = liver.merge(
        generic,
        on=["fold_index", "true_kinase_node_id", "site_node_id"],
        suffixes=("_liver", "_generic"),
    )

    liver_ttk_cases = merged[merged["top1_predicted_kinase_liver"] == TTK].copy()

    n_same_rows_generic_ttk = (liver_ttk_cases["top1_predicted_kinase_generic"] == TTK).sum()
    n_same_rows_generic_not_ttk = (liver_ttk_cases["top1_predicted_kinase_generic"] != TTK).sum()

    # Optional: among true TTK substrates, how often does liver/generic rank TTK top1
    ttk_true_merged = merged[merged["true_kinase_node_id"] == TTK].copy()
    n_true_ttk_liver_top1 = (ttk_true_merged["top1_predicted_kinase_liver"] == TTK).sum()
    n_true_ttk_generic_top1 = (ttk_true_merged["top1_predicted_kinase_generic"] == TTK).sum()

    print("TTK summary")
    print("=" * 60)
    print(f"Known TTK substrate rows in results: {n_ttk_true_rows}")
    print(f"Known TTK substrate sites in results: {n_ttk_true_sites}")
    print()
    print(f"Liver predictions with TTK as top1 (rows): {n_liver_ttk_top1_rows}")
    print(f"Liver predictions with TTK as top1 (unique sites): {n_liver_ttk_top1_sites}")
    print()
    print("Among rows where liver predicts TTK as top1:")
    print(f"Generic also predicts TTK as top1: {n_same_rows_generic_ttk}")
    print(f"Generic predicts a different top1: {n_same_rows_generic_not_ttk}")
    print()
    print("Among true TTK substrate rows:")
    print(f"Liver predicts TTK as top1: {n_true_ttk_liver_top1}")
    print(f"Generic predicts TTK as top1: {n_true_ttk_generic_top1}")
    print()

    if not liver_ttk_cases.empty:
        out_cols = [
            "fold_index",
            "true_kinase_node_id",
            "site_node_id",
            "top1_predicted_kinase_liver",
            "top1_score_liver",
            "held_out_rank_liver",
            "best_true_rank_liver",
            "top1_predicted_kinase_generic",
            "top1_score_generic",
            "held_out_rank_generic",
            "best_true_rank_generic",
        ]
        available = [c for c in out_cols if c in liver_ttk_cases.columns]
        preview = liver_ttk_cases[available].sort_values(
            by=[c for c in ["site_node_id", "fold_index"] if c in liver_ttk_cases.columns]
        )
        print("Preview of cases where liver predicts TTK as top1:")
        print(preview.head(20).to_string(index=False))

        preview.to_csv(OUTPUT_DIR / "ttk_focus_cases.csv", index=False)
        print()
        print(f"Wrote: {OUTPUT_DIR / 'ttk_focus_cases.csv'}")


if __name__ == "__main__":
    main()