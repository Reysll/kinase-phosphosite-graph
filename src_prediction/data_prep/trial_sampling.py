from __future__ import annotations

import numpy as np
import pandas as pd


def sample_trials(
    positive_edges: pd.DataFrame,
    max_trials: int,
    random_state: int = 42,
) -> pd.DataFrame:
    trials = positive_edges.copy().reset_index(drop=True)

    if max_trials < len(trials):
        trials = trials.sample(n=max_trials, random_state=random_state).reset_index(drop=True)

    trials = trials.reset_index(drop=True)
    trials["trial_index"] = np.arange(1, len(trials) + 1)
    return trials
