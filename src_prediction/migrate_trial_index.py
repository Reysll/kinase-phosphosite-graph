"""
One-time migration: rename fold_index -> trial_index in existing result and
frozen-trial CSV files that were generated before the terminology change.
Safe to re-run: skips files that already have trial_index.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src_prediction.config import PRED_OUTPUTS_DIR


FILES = [
    PRED_OUTPUTS_DIR / "generic_multi_kinase" / "results.csv.gz",
    PRED_OUTPUTS_DIR / "liver_multi_kinase" / "results.csv.gz",
    PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz",
    PRED_OUTPUTS_DIR / "frozen_trials_10.csv.gz",
    PRED_OUTPUTS_DIR / "frozen_trials_50.csv.gz",
]


def migrate_file(path: Path) -> None:
    if not path.exists():
        print(f"  skip (not found): {path.name}")
        return

    df = pd.read_csv(path)

    if "trial_index" in df.columns:
        print(f"  already migrated:  {path.name}")
        return

    if "fold_index" not in df.columns:
        print(f"  no fold_index:     {path.name}")
        return

    df = df.rename(columns={"fold_index": "trial_index"})
    df.to_csv(path, index=False, compression="gzip")
    print(f"  migrated:          {path.name}  ({len(df):,} rows)")


def main() -> None:
    print("Migrating fold_index -> trial_index in existing result files...\n")
    for f in FILES:
        migrate_file(f)
    print("\nDone.")


if __name__ == "__main__":
    main()
