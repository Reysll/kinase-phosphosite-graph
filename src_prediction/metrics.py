from __future__ import annotations

from typing import Iterable

import pandas as pd


def top_k_accuracy(results_df: pd.DataFrame, k: int) -> float:
    valid = results_df["rank_of_true_kinase"].dropna()
    if len(valid) == 0:
        return 0.0
    return float((valid <= k).mean())


def summarize_results(results_df: pd.DataFrame, ks: Iterable[int] = (1, 5, 10, 20)) -> pd.DataFrame:
    rows = []

    n_folds = len(results_df)
    n_ranked = results_df["rank_of_true_kinase"].notna().sum()
    n_missing_rank = results_df["rank_of_true_kinase"].isna().sum()

    for k in ks:
        rows.append(
            {
                "metric": f"top_{k}_accuracy",
                "value": top_k_accuracy(results_df, k),
            }
        )

    rows.extend(
        [
            {"metric": "n_folds", "value": float(n_folds)},
            {"metric": "n_ranked", "value": float(n_ranked)},
            {"metric": "n_missing_rank", "value": float(n_missing_rank)},
            {
                "metric": "mean_rank",
                "value": float(results_df["rank_of_true_kinase"].dropna().mean()) if n_ranked > 0 else float("nan"),
            },
            {
                "metric": "median_rank",
                "value": float(results_df["rank_of_true_kinase"].dropna().median()) if n_ranked > 0 else float("nan"),
            },
        ]
    )

    return pd.DataFrame(rows)


def summarize_results_text(results_df: pd.DataFrame, ks: Iterable[int] = (1, 5, 10, 20)) -> str:
    summary_df = summarize_results(results_df, ks=ks)

    lines = ["Baseline node2vec + cosine similarity results", "=========================================="]
    for row in summary_df.itertuples(index=False):
        lines.append(f"{row.metric}: {row.value}")

    return "\n".join(lines)