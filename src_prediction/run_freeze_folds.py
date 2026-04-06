from __future__ import annotations

from pathlib import Path

from src_prediction.config import LIVER_POSITIVE_EDGES_OUT, PRED_OUTPUTS_DIR
from src_prediction.fold_sampling import sample_folds
from src_prediction.io_utils import read_csv_auto, write_csv_gz


def main() -> None:
    positive_edges = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)

    max_folds = 50
    random_state = 42

    frozen = sample_folds(
        positive_edges=positive_edges,
        max_folds=max_folds,
        random_state=random_state,
    )

    out_path = PRED_OUTPUTS_DIR / "frozen_folds_50.csv.gz"
    write_csv_gz(frozen, out_path)

    print(f"Frozen folds: {len(frozen):,}")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()