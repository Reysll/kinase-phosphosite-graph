# src_prediction/make_multi_kinase_fold_set_10.py

from pathlib import Path
import pandas as pd


INPUT_PATH = Path("outputs_prediction/frozen_folds_multi_kinase_sites.csv.gz")
OUTPUT_PATH = Path("outputs_prediction/frozen_folds_multi_kinase_sites_10.csv.gz")


def main() -> None:
    if not INPUT_PATH.exists():
        raise FileNotFoundError(
            f"Missing input file: {INPUT_PATH}\n"
            f"Run python -m src_prediction.make_multi_kinase_fold_set first."
        )

    print(f"Reading: {INPUT_PATH}")
    df = pd.read_csv(INPUT_PATH)

    if "site_node_id" not in df.columns:
        raise ValueError(
            "Expected column 'site_node_id' was not found in the input file."
        )

    # Keep one row per site so the 10 folds cover 10 different sites
    df10 = df.drop_duplicates(subset=["site_node_id"]).head(10).copy()

    # Renumber fold_index from 1 to 10
    if "fold_index" in df10.columns:
        df10["fold_index"] = range(1, len(df10) + 1)

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    df10.to_csv(OUTPUT_PATH, index=False, compression="gzip")

    print()
    print(f"Rows written: {len(df10)}")
    print(f"Unique sites written: {df10['site_node_id'].nunique()}")
    print(f"Wrote: {OUTPUT_PATH}")
    print()
    print("Preview:")
    print(df10.to_string(index=False))


if __name__ == "__main__":
    main()