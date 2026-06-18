"""
Kinase-Substrate Enrichment Analysis (KSEA).

Formula from Wiredja et al. 2017 (Bioinformatics 33:3489):
  score_k = (s_bar - p_bar) / (m * delta)

where:
  s_bar  = mean log2(FC) of kinase k's substrates found in the dataset
  p_bar  = mean log2(FC) of all phosphosites in the dataset
  m      = number of kinase k's substrates found in the dataset
  delta  = standard deviation of log2(FC) across all phosphosites

FC = log2(mean_tumor_abundance / mean_ctrl_abundance) per phosphosite.
Positive score -> substrates are hyper-phosphorylated in tumor vs control.

P-values: one-tailed probability under N(0,1) that score is more extreme,
followed by Benjamini-Hochberg FDR correction (Wiredja 2017, Section 2.1).
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Dict, Set

import numpy as np
import pandas as pd
from scipy import stats


_PHOSPHO_SITE_RE = re.compile(r"([STY])\s*([0-9]+)")


def load_site_fc_values(excel_path: Path, min_identified: int = 6) -> pd.Series:
    """
    Parse the Phosphorylation sheet from the liver proteomics Excel file.

    Returns a Series indexed by site_node_id (SITE:GENE-S123 format) with
    log2(mean_tumor / mean_ctrl) as the value. Sites with insufficient
    measurements or zero/NaN control mean are excluded.
    """
    df = pd.read_excel(excel_path, sheet_name="Phosphorylation", header=1)

    ctrl_cols = [c for c in df.columns if "Control" in str(c)]
    samp_cols = [c for c in df.columns if "Sample" in str(c)]

    if not ctrl_cols or not samp_cols:
        raise ValueError("Could not find Control/Sample columns in Phosphorylation sheet.")

    ctrl_mat = df[ctrl_cols].apply(pd.to_numeric, errors="coerce")
    samp_mat = df[samp_cols].apply(pd.to_numeric, errors="coerce")

    # Build site node IDs matching the graph format: SITE:GENE-S123
    gene_list = [str(g) for g in df["Gene Symbol"].tolist()]
    mod_list = df["Modifications in Master Proteins"].tolist()
    site_ids: list[str] = []
    for gene, mod in zip(gene_list, mod_list):
        gene = gene.strip()
        hits = _PHOSPHO_SITE_RE.findall(str(mod)) if isinstance(mod, str) else []
        if len(hits) == 1 and gene and gene != "nan":
            aa, pos = hits[0]
            site_ids.append(f"SITE:{gene}-{aa}{pos}")
        else:
            site_ids.append("")

    valid = pd.Series([s != "" for s in site_ids], index=df.index)
    ctrl_ok = ctrl_mat.notna().sum(axis=1) >= 1
    samp_ok = samp_mat.notna().sum(axis=1) >= min_identified
    keep = valid & ctrl_ok & samp_ok

    ctrl_mean = ctrl_mat.loc[keep].mean(axis=1, skipna=True).replace(0, np.nan)
    samp_mean = samp_mat.loc[keep].mean(axis=1, skipna=True)

    log2_fc = np.log2(samp_mean / ctrl_mean)
    log2_fc.index = pd.Index([site_ids[i] for i in keep[keep].index])

    log2_fc = log2_fc.dropna()
    # Collapse duplicate site IDs by averaging
    log2_fc = log2_fc.groupby(log2_fc.index).mean()

    return log2_fc


def compute_ksea_scores(
    fc_values: pd.Series,
    substrate_sets: Dict[str, Set[str]],
    min_substrates: int = 3,
) -> pd.DataFrame:
    """
    Compute KSEA scores for all kinases using the Wiredja 2017 formula.

    score_k = (s_bar - p_bar) / (m * delta)

    where delta is global SD and m is the substrate count — NOT delta/sqrt(m).
    This penalizes kinases with larger substrate sets, making the score
    sensitive to the per-substrate FC enrichment rather than coverage.

    Args:
        fc_values: log2FC per site (indexed by site_node_id).
        substrate_sets: kinase_node_id -> set of site_node_ids.
        min_substrates: minimum substrates present in fc_values to score a kinase.

    Returns DataFrame with columns:
        kinase_node_id, kinase, n_substrates_in_dataset,
        mean_substrate_fc, global_mean_fc, ksea_score, p_value, adj_p_value
    Sorted by ksea_score descending.
    """
    global_mean = fc_values.mean()
    global_sd = fc_values.std()

    rows = []
    for kinase_id, substrates in substrate_sets.items():
        present = [s for s in substrates if s in fc_values.index]
        n = len(present)
        if n < min_substrates:
            continue

        mean_sub = fc_values[present].mean()
        # Wiredja 2017: denominator is m * delta (NOT delta/sqrt(m))
        score = (mean_sub - global_mean) / (n * global_sd)

        # One-tailed p-value: P(Z >= |score|) for positive, P(Z <= score) for negative
        p_value = stats.norm.sf(abs(score))

        rows.append(
            {
                "kinase_node_id": kinase_id,
                "kinase": kinase_id.split(":", 1)[-1],
                "n_substrates_in_dataset": n,
                "mean_substrate_fc": mean_sub,
                "global_mean_fc": global_mean,
                "ksea_score": score,
                "p_value": p_value,
            }
        )

    if not rows:
        return pd.DataFrame(
            columns=[
                "kinase_node_id",
                "kinase",
                "n_substrates_in_dataset",
                "mean_substrate_fc",
                "global_mean_fc",
                "ksea_score",
                "p_value",
                "adj_p_value",
            ]
        )

    result = pd.DataFrame(rows).sort_values("ksea_score", ascending=False).reset_index(drop=True)

    # Benjamini-Hochberg FDR correction
    if hasattr(stats, "false_discovery_control"):
        adj_p = stats.false_discovery_control(result["p_value"].values, method="bh")
    else:
        adj_p = _bh_correction(result["p_value"].values)
    result["adj_p_value"] = adj_p

    return result


def _bh_correction(p_values: np.ndarray) -> tuple:
    """Benjamini-Hochberg FDR correction for scipy versions without false_discovery_control."""
    n = len(p_values)
    order = np.argsort(p_values)
    ranked_p = p_values[order]
    adj = np.minimum(1.0, ranked_p * n / (np.arange(1, n + 1)))
    # Enforce monotonicity (cumulative min from right)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    result = np.empty(n)
    result[order] = adj
    return None, result, None, None
