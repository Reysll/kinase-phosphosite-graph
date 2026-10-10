"""
G2 -- KSEA p-value convention check (two-tailed recomputation).

READ-ONLY. Reads outputs_prediction/ksea/ksea_psp.csv and writes only into
thesis/analysis/out/. It does not modify any project code, data or results.

Issue
-----
src_prediction/analysis/ksea.py computes

    p = scipy.stats.norm.sf(abs(score))

which is a ONE-tailed tail probability evaluated at the ABSOLUTE score. Because
the sign of the score is discarded, the test is effectively two-sided in intent
but omits the factor of 2, so every p-value is half of the proper two-tailed
value. The published count of 17/64 significant kinases rests on that
convention.

This script recomputes

    p_two_tailed = 2 * norm.sf(abs(score))

and applies the SAME multiple-testing correction as the original
(Benjamini-Hochberg via scipy.stats.false_discovery_control, method="bh"),
then reports both counts side by side.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.ksea_twotailed

Seed: none required (deterministic).
"""

from __future__ import annotations

import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

KSEA_CSV = Path("outputs_prediction/ksea/ksea_psp.csv")
OUT = Path("thesis/analysis/out")
ALPHA = 0.05


def git_hash_for(path: Path) -> str:
    try:
        out = subprocess.run(
            ["git", "log", "-1", "--format=%H", "--", str(path)],
            capture_output=True, text=True, check=True,
        ).stdout.strip()
        return out or "UNTRACKED/UNCOMMITTED"
    except Exception:
        return "UNKNOWN"


def bh(p: np.ndarray) -> np.ndarray:
    return np.asarray(stats.false_discovery_control(p, method="bh"))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(KSEA_CSV)

    z = df["ksea_score"].to_numpy(dtype=float)
    p_one = stats.norm.sf(np.abs(z))          # as implemented in ksea.py
    p_two = 2.0 * stats.norm.sf(np.abs(z))    # proper two-tailed

    df["p_recomputed_one_tailed"] = p_one
    df["adj_p_recomputed_one_tailed"] = bh(p_one)
    df["p_two_tailed"] = p_two
    df["adj_p_two_tailed"] = bh(p_two)
    df["direction"] = np.where(z > 0, "up in tumor", "down in tumor")

    n_total = len(df)
    n_orig = int((df["adj_p_value"] < ALPHA).sum())
    n_one = int((df["adj_p_recomputed_one_tailed"] < ALPHA).sum())
    n_two = int((df["adj_p_two_tailed"] < ALPHA).sum())

    sig_one = set(df.loc[df["adj_p_recomputed_one_tailed"] < ALPHA, "kinase"])
    sig_two = set(df.loc[df["adj_p_two_tailed"] < ALPHA, "kinase"])
    lost = sorted(sig_one - sig_two)

    two = df[df["adj_p_two_tailed"] < ALPHA].sort_values("ksea_score", ascending=False)

    lines: list[str] = []
    add = lines.append
    add("G2 -- KSEA p-value convention: one-tailed-on-|z| vs two-tailed")
    add("=" * 78)
    add("")
    add(f"Source file   : {str(KSEA_CSV).replace(chr(92), '/')}")
    add(f"git commit    : {git_hash_for(KSEA_CSV)}")
    add(f"Kinases scored: {n_total}   (>= 3 substrates in the FC dataset)")
    add(f"Correction    : Benjamini-Hochberg, alpha = {ALPHA}")
    add("")
    add("Significant kinase counts")
    add("-" * 78)
    add(f"  stored adj_p_value (as published)          : {n_orig} / {n_total}")
    add(f"  recomputed p = norm.sf(|z|)   [reproduces] : {n_one} / {n_total}")
    add(f"  recomputed p = 2*norm.sf(|z|) [two-tailed] : {n_two} / {n_total}")
    add("")
    add(f"Kinases that lose significance under the two-tailed convention ({len(lost)}):")
    add("  " + (", ".join(lost) if lost else "none"))
    add("")
    add(f"Two-tailed significant kinases ({n_two})")
    add("-" * 78)
    show = two[["kinase", "n_substrates_in_dataset", "mean_substrate_fc",
                "ksea_score", "p_two_tailed", "adj_p_two_tailed", "direction"]].copy()
    show["mean_substrate_fc"] = show["mean_substrate_fc"].map("{:.3f}".format)
    show["ksea_score"] = show["ksea_score"].map("{:+.4f}".format)
    for c in ("p_two_tailed", "adj_p_two_tailed"):
        show[c] = show[c].map(lambda v: f"{v:.3e}" if v < 1e-3 else f"{v:.4f}")
    add(show.to_string(index=False))
    add("")
    add("Direction split (two-tailed significant): "
        f"up = {int((two['ksea_score'] > 0).sum())}, "
        f"down = {int((two['ksea_score'] < 0).sum())}")

    text = "\n".join(lines)
    (OUT / "ksea_twotailed.txt").write_text(text, encoding="utf-8")
    df.sort_values("ksea_score", ascending=False).to_csv(
        OUT / "ksea_twotailed.csv", index=False
    )
    print(text)


if __name__ == "__main__":
    main()
