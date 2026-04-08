from __future__ import annotations

from typing import Iterable

import pandas as pd


def top_k_accuracy_from_col(results_df: pd.DataFrame, rank_col: str, k: int) -> float:
    valid = results_df[rank_col].dropna()
    if len(valid) == 0:
        return 0.0
    return float((valid <= k).mean())


def summarize_results(
    results_df: pd.DataFrame,
    ks: Iterable[int] = (1, 5, 10, 20),
) -> pd.DataFrame:
    rows = []

    n_folds = len(results_df)

    for rank_col, prefix in [
        ("held_out_kinase_rank", "held_out"),
        ("best_true_kinase_rank", "best_true"),
    ]:
        n_ranked = results_df[rank_col].notna().sum()
        n_missing_rank = results_df[rank_col].isna().sum()

        for k in ks:
            rows.append(
                {
                    "metric": f"{prefix}_top_{k}_accuracy",
                    "value": top_k_accuracy_from_col(results_df, rank_col, k),
                }
            )

        rows.extend(
            [
                {"metric": f"{prefix}_n_ranked", "value": float(n_ranked)},
                {"metric": f"{prefix}_n_missing_rank", "value": float(n_missing_rank)},
                {
                    "metric": f"{prefix}_mean_rank",
                    "value": float(results_df[rank_col].dropna().mean()) if n_ranked > 0 else float("nan"),
                },
                {
                    "metric": f"{prefix}_median_rank",
                    "value": float(results_df[rank_col].dropna().median()) if n_ranked > 0 else float("nan"),
                },
            ]
        )

    rows.append({"metric": "n_folds", "value": float(n_folds)})
    rows.append(
        {
            "metric": "mean_n_true_kinases_per_site",
            "value": float(results_df["n_true_kinases_for_site"].mean()) if n_folds > 0 else float("nan"),
        }
    )

    return pd.DataFrame(rows)


def summarize_results_text(
    results_df: pd.DataFrame,
    ks: Iterable[int] = (1, 5, 10, 20),
) -> str:
    summary_df = summarize_results(results_df, ks=ks)

    lines = [
        "Baseline node2vec results with multi-true-kinase evaluation",
        "===========================================================",
    ]
    for row in summary_df.itertuples(index=False):
        lines.append(f"{row.metric}: {row.value}")

    return "\n".join(lines)