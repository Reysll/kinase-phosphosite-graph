from __future__ import annotations

import numpy as np
import pandas as pd


def sample_folds(
    positive_edges: pd.DataFrame,
    max_folds: int,
    random_state: int = 42,
) -> pd.DataFrame:
    folds = positive_edges.copy().reset_index(drop=True)

    if max_folds < len(folds):
        folds = folds.sample(n=max_folds, random_state=random_state).reset_index(drop=True)

    folds = folds.reset_index(drop=True)
    folds["fold_index"] = np.arange(1, len(folds) + 1)
    return folds