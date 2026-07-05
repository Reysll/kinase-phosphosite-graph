# src_prediction/check_ttk_from_existing_debug_results.py

from pathlib import Path
import pandas as pd

OUTPUT_DIR = Path("outputs_prediction")
COMPARE_FILE = OUTPUT_DIR / "comparison_generic_vs_liver_fc_corr.csv.gz"
TTK = "PROTEIN:TTK"


def main():
    if not COMPARE_FILE.exists():
        raise FileNotFoundError(
            f"Missing file: {COMPARE_FILE}\n"
            f"This script uses your existing generic vs liver debug comparison output."
        )

    df = pd.read_csv(COMPARE_FILE)

    # In how many liver-specific predictions is TTK top1
    liver_ttk = df[df["liver_top1"] == TTK].copy()

    # Among those, was generic also top1 = TTK
    generic_also_ttk = liver_ttk[liver_ttk["generic_top1"] == TTK].copy()
    generic_not_ttk = liver_ttk[liver_ttk["generic_top1"] != TTK].copy()

    print("TTK in existing generic vs liver debug comparison")
    print("=" * 60)
    print(f"Total rows in comparison: {len(df)}")
    print(f"Liver predicts TTK as top1 in: {len(liver_ttk)} row(s)")
    print(f"Generic also predicts TTK as top1 in: {len(generic_also_ttk)} of those row(s)")
    print(f"Generic predicts a different top1 in: {len(generic_not_ttk)} of those row(s)")
    print()

    if not liver_ttk.empty:
        print("Rows where liver predicts TTK as top1:")
        cols = [
            "fold_index",
            "true_kinase_node_id",
            "site_node_id",
            "generic_top1",
            "generic_top1_score",
            "generic_held_out_rank",
            "liver_top1",
            "liver_top1_score",
            "liver_held_out_rank",
            "held_out_rank_delta",
            "best_true_rank_delta",
        ]
        cols = [c for c in cols if c in liver_ttk.columns]
        print(liver_ttk[cols].to_string(index=False))

        out_path = OUTPUT_DIR / "ttk_from_existing_debug_results.csv"
        liver_ttk[cols].to_csv(out_path, index=False)
        print()
        print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()