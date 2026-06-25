from __future__ import annotations

from src_prediction.core.config import LIVER_POSITIVE_EDGES_OUT, PRED_OUTPUTS_DIR
from src_prediction.data_prep.trial_sampling import sample_trials
from src_prediction.core.io_utils import read_csv_auto, write_csv_gz


def main() -> None:
    positive_edges = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)

    max_trials = 10
    random_state = 42

    frozen = sample_trials(
        positive_edges=positive_edges,
        max_trials=max_trials,
        random_state=random_state,
    )

    out_path = PRED_OUTPUTS_DIR / "frozen_trials_10.csv.gz"
    write_csv_gz(frozen, out_path)

    print(f"Frozen LOO trials: {len(frozen):,}")
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
