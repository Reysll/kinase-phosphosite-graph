"""
Result 2 (Dr. Ayati, agreed 2026-06-24): per-kinase rank-distribution
KS-test between two networks, across the 3,228 FC-valid sites.

Two network pairs, per Salvador's 2026-06-24 instruction (run both, see
what Dr. Ayati thinks):
  1. control vs liver FC   -- Dr. Ayati's original request (healthy vs
     fold-change-disease contrast).
  2. cancer vs control     -- the alternative "cleaner" healthy-vs-disease
     contrast proposed in email_reply.txt (raw tumor abundance vs raw
     healthy abundance, no FC normalization mixed in).

Run across all three embeddings (node2vec, spectral, concat) since the
inference-pass output already exists for all of them.

Neither pair involves "generic", so both sides of every comparison have full
3,228-site coverage (no coverage-asymmetry caveat here, unlike Result 1).

Consumes outputs_prediction/inference_all/{network}_{embedding}/ranks.csv.gz.

IMPORTANT CAVEAT (found 2026-06-25, while building this script): for most
kinases, rank_in_site is nearly CONSTANT across all 3,228 sites within one
network (e.g. CDK20 is rank 349 at literally every control site, rank 350 at
literally every cancer site -- zero within-network variance). This is the
global model's hub bias showing up in a new way: for non-true-kinase pairs,
the kinase's own embedding dominates its score far more than the site's
embedding does, so its rank barely depends on which site is being scored.
A KS-test on two near-degenerate (point-mass-like) distributions that differ
by even one rank position gives ks_statistic ~= 1.0 trivially -- this is why
~414/415 kinases come back "significant" below. That number is statistically
correct but NOT a meaningful biological finding by itself.
Fix: every kinase's within-network rank std-dev is reported
(std_rank_{a,b}), and an `informative` flag marks kinases where at least one
side has std_rank > MIN_INFORMATIVE_STD (i.e. its rank genuinely varies by
site, not just a fixed per-kinase baseline). The summary's top-significant
list is restricted to informative kinases; the per-pair CSVs keep everything
so the degenerate cases are visible, not hidden.

Outputs (outputs_prediction/analysis/):
  ks_test_{a}_vs_{b}_{embedding}.csv   -- per-kinase KS-test results (all kinases)
  ks_test_summary.txt                   -- informative + significant kinases, all pairs/embeddings
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from src_prediction.core.config import PRED_OUTPUTS_DIR
from src_prediction.core.io_utils import write_text

INFERENCE_ALL_DIR = PRED_OUTPUTS_DIR / "inference_all"
ANALYSIS_DIR = PRED_OUTPUTS_DIR / "analysis"

PAIRS = [("control", "liver_fc"), ("cancer", "control")]
EMBEDDINGS = ["node2vec", "spectral", "concat"]
MIN_N_PER_SIDE = 30  # need enough sites per kinase for a meaningful KS-test
# A kinase's rank is "informative" if it actually varies by site in at least
# one of the two networks being compared (std_rank > this threshold). Below
# this, the rank is effectively a fixed per-kinase baseline -- see module
# docstring caveat.
MIN_INFORMATIVE_STD = 1.0


def _clean_kinase(node_id: str) -> str:
    return str(node_id).split(":", 1)[-1]


def _bh_correction(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR correction (fallback for older scipy)."""
    n = len(p_values)
    order = np.argsort(p_values)
    ranked_p = p_values[order]
    adj = np.minimum(1.0, ranked_p * n / (np.arange(1, n + 1)))
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    result = np.empty(n)
    result[order] = adj
    return result


def _load_ranks(network: str, embedding: str) -> pd.DataFrame | None:
    path = INFERENCE_ALL_DIR / f"{network}_{embedding}" / "ranks.csv.gz"
    if not path.exists():
        return None
    return pd.read_csv(path)[["kinase_node_id", "rank_in_site"]]


def run_ks_test(network_a: str, network_b: str, embedding: str) -> pd.DataFrame:
    df_a = _load_ranks(network_a, embedding)
    df_b = _load_ranks(network_b, embedding)
    if df_a is None or df_b is None:
        return pd.DataFrame()

    groups_a = df_a.groupby("kinase_node_id")["rank_in_site"].apply(np.array)
    groups_b = df_b.groupby("kinase_node_id")["rank_in_site"].apply(np.array)
    common_kinases = sorted(set(groups_a.index) & set(groups_b.index))

    rows = []
    for kinase in common_kinases:
        ranks_a = groups_a[kinase]
        ranks_b = groups_b[kinase]
        if len(ranks_a) < MIN_N_PER_SIDE or len(ranks_b) < MIN_N_PER_SIDE:
            continue

        ks_stat, p_value = stats.ks_2samp(ranks_a, ranks_b)
        std_a, std_b = ranks_a.std(), ranks_b.std()
        rows.append(
            {
                "kinase_node_id": kinase,
                "kinase": _clean_kinase(kinase),
                f"n_{network_a}": len(ranks_a),
                f"n_{network_b}": len(ranks_b),
                f"mean_rank_{network_a}": ranks_a.mean(),
                f"mean_rank_{network_b}": ranks_b.mean(),
                f"median_rank_{network_a}": np.median(ranks_a),
                f"median_rank_{network_b}": np.median(ranks_b),
                f"std_rank_{network_a}": std_a,
                f"std_rank_{network_b}": std_b,
                "ks_statistic": ks_stat,
                "p_value": p_value,
                "shift": (
                    f"better in {network_b}"
                    if np.median(ranks_b) < np.median(ranks_a)
                    else f"better in {network_a}"
                ),
                "informative": bool(std_a > MIN_INFORMATIVE_STD or std_b > MIN_INFORMATIVE_STD),
            }
        )

    if not rows:
        return pd.DataFrame()

    result = pd.DataFrame(rows)
    if hasattr(stats, "false_discovery_control"):
        result["adj_p_value"] = stats.false_discovery_control(result["p_value"].values, method="bh")
    else:
        result["adj_p_value"] = _bh_correction(result["p_value"].values)

    return result.sort_values(["adj_p_value", "ks_statistic"], ascending=[True, False]).reset_index(
        drop=True
    )


def main() -> None:
    ANALYSIS_DIR.mkdir(parents=True, exist_ok=True)

    summary_lines = [
        "Result 2 -- Per-Kinase Rank-Distribution KS-Test",
        "=" * 70,
        "",
        "scipy.stats.ks_2samp on each kinase's rank_in_site distribution across",
        "the 3,228 FC-valid sites, comparing two networks. BH-FDR corrected",
        f"across kinases. Kinases need >= {MIN_N_PER_SIDE} scored sites per side.",
        "",
        "CAVEAT: most kinases' rank_in_site is nearly constant across all sites",
        "within one network (hub bias -- the kinase's own embedding dominates",
        "its score far more than the site's does), so a 1-rank shift between",
        "networks gives ks_statistic ~= 1.0 trivially. 'Significant' below is",
        f"restricted to kinases flagged `informative` (std_rank > {MIN_INFORMATIVE_STD}",
        "in at least one network -- i.e. its rank genuinely varies by site).",
        "Full per-pair CSVs include ALL kinases (informative and not) with",
        "std_rank columns so the degenerate cases stay visible.",
        "",
    ]

    for network_a, network_b in PAIRS:
        for embedding in EMBEDDINGS:
            result = run_ks_test(network_a, network_b, embedding)
            tag = f"{network_a}_vs_{network_b}_{embedding}"
            if result.empty:
                summary_lines.append(f"{tag}: no data, skipped.\n")
                continue

            csv_path = ANALYSIS_DIR / f"ks_test_{tag}.csv"
            result.to_csv(csv_path, index=False)
            print(f"Wrote: {csv_path}")

            n_sig = int((result["adj_p_value"] < 0.05).sum())
            n_informative = int(result["informative"].sum())
            n_sig_informative = int(((result["adj_p_value"] < 0.05) & result["informative"]).sum())
            summary_lines += [
                f"--- {network_a} vs {network_b} -- {embedding} ---",
                f"Kinases tested: {len(result)}   "
                f"Significant (adj_p < 0.05): {n_sig}   "
                f"Informative (rank varies by site): {n_informative}   "
                f"Significant AND informative: {n_sig_informative}",
            ]
            top_sig = result[(result["adj_p_value"] < 0.05) & result["informative"]].head(15)
            if top_sig.empty:
                summary_lines.append("  No significant + informative kinases.")
            else:
                for _, row in top_sig.iterrows():
                    summary_lines.append(
                        f"  {row['kinase']:<14}  ks={row['ks_statistic']:.3f}  "
                        f"adj_p={row['adj_p_value']:.2e}  "
                        f"median {network_a}={row[f'median_rank_{network_a}']:.0f} "
                        f"(std={row[f'std_rank_{network_a}']:.1f})  "
                        f"median {network_b}={row[f'median_rank_{network_b}']:.0f} "
                        f"(std={row[f'std_rank_{network_b}']:.1f})  "
                        f"({row['shift']})"
                    )
            summary_lines.append("")

    write_text("\n".join(summary_lines), ANALYSIS_DIR / "ks_test_summary.txt")
    print(f"Wrote: {ANALYSIS_DIR / 'ks_test_summary.txt'}")


if __name__ == "__main__":
    main()
