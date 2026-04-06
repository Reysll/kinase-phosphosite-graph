from __future__ import annotations

from pathlib import Path

import pandas as pd

from src_prediction.io_utils import write_csv_gz, write_text


def experiment_paths(base_dir: Path, experiment_name: str) -> dict[str, Path]:
    exp_dir = base_dir / experiment_name
    exp_dir.mkdir(parents=True, exist_ok=True)

    return {
        "dir": exp_dir,
        "results": exp_dir / "results.csv.gz",
        "metrics": exp_dir / "metrics.csv.gz",
        "summary": exp_dir / "summary.txt",
    }


def write_experiment_outputs(
    experiment_name: str,
    base_dir: Path,
    results_df: pd.DataFrame,
    metrics_df: pd.DataFrame,
    summary_text: str,
) -> dict[str, Path]:
    paths = experiment_paths(base_dir, experiment_name)

    write_csv_gz(results_df, paths["results"])
    write_csv_gz(metrics_df, paths["metrics"])
    write_text(summary_text, paths["summary"])

    return paths