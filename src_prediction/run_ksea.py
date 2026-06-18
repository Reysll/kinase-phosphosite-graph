"""
Run KSEA comparing kinase activity under three substrate set sources:
  1. PSP known KSAs (from generic graph phosphorylates edges)
  2. Top-1 model predictions — generic spectral (which sites does generic predict each kinase for?)
  3. Top-1 model predictions — liver FC spectral

The top-1 prediction approach: for each trial, top1_predicted_kinase is the globally
top-scoring kinase. Collecting all (kinase, site) pairs where a kinase was predicted
as top-1 gives its "predicted substrate set" from the model's perspective.

Outputs written to outputs_prediction/ksea/
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src_prediction.config import GENERIC_EDGES, PRED_OUTPUTS_DIR
from src_prediction.io_utils import read_csv_auto, write_text
from src_prediction.ksea import compute_ksea_scores, load_site_fc_values

LIVER_EXCEL = Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")
KSEA_DIR = PRED_OUTPUTS_DIR / "ksea"
MIN_SUBSTRATES = 3


def _build_psp_substrate_sets(generic_edges: pd.DataFrame) -> dict[str, set[str]]:
    psp = generic_edges[generic_edges["relation"] == "phosphorylates"]
    result: dict[str, set[str]] = {}
    for kinase, group in psp.groupby("source"):
        result[kinase] = set(group["target"].tolist())
    return result


def _build_top1_substrate_sets(results_df: pd.DataFrame) -> dict[str, set[str]]:
    """
    For each kinase, collect all sites where the model predicted it as top-1.
    Uses the top1_predicted_kinase column from LOO results.
    """
    if "top1_predicted_kinase" not in results_df.columns:
        return {}
    result: dict[str, set[str]] = {}
    for _, row in results_df[["site_node_id", "top1_predicted_kinase"]].dropna().iterrows():
        kinase = row["top1_predicted_kinase"]
        site = row["site_node_id"]
        result.setdefault(kinase, set()).add(site)
    return result


def _fmt_table(df: pd.DataFrame, title: str, top_n: int = 20) -> list[str]:
    out = [title, "-" * len(title)]
    if df.empty:
        out += ["  — no data", ""]
        return out
    out.append(
        f"  {'Kinase':<20}  {'n_subs':>6}  {'mean_FC':>8}  {'score':>10}  {'p_val':>8}  {'adj_p':>8}"
    )
    out.append("  " + "-" * 68)
    for _, row in df.head(top_n).iterrows():
        p = row.get("p_value", float("nan"))
        ap = row.get("adj_p_value", float("nan"))
        out.append(
            f"  {row['kinase']:<20}  {int(row['n_substrates_in_dataset']):>6}  "
            f"{row['mean_substrate_fc']:>8.3f}  {row['ksea_score']:>10.4f}  "
            f"{p:>8.4f}  {ap:>8.4f}"
        )
    if len(df) > top_n:
        out.append(f"  ... ({len(df) - top_n} more kinases)")
    out.append("")
    return out


def _build_report(
    fc_values: pd.Series,
    psp_scores: pd.DataFrame,
    generic_scores: pd.DataFrame,
    liver_scores: pd.DataFrame,
) -> str:
    lines = [
        "KSEA Comparison Report — Generic vs Liver FC Spectral",
        "=" * 70,
        "",
        "Formula (Wiredja et al. 2017, Bioinformatics 33:3489):",
        "  score_k = (s_bar - p_bar) / (m * delta)",
        "    s_bar  = mean log2(FC) of kinase substrates in dataset",
        "    p_bar  = mean log2(FC) of all phosphosites (global mean)",
        "    m      = number of kinase substrates found in dataset",
        "    delta  = standard deviation of all log2(FC) values",
        "  FC = log2(mean_tumor / mean_ctrl) per phosphosite.",
        "  Positive score -> substrates hyper-phosphorylated in tumor.",
        "  NOTE: m in the denominator penalizes kinases with many substrates.",
        "  Kinases with few high-FC substrates score higher than kinases with",
        "  many substrates at moderate FC, even if the latter is more significant.",
        "  P-values: one-tailed N(0,1), Benjamini-Hochberg FDR correction.",
        "",
        f"FC dataset: {len(fc_values):,} phosphosites",
        f"Global mean log2FC: {fc_values.mean():.3f}   SD: {fc_values.std():.3f}",
        f"Min substrates to score a kinase: {MIN_SUBSTRATES}",
        "",
    ]

    lines += _fmt_table(psp_scores, "Top kinases — PSP known substrates")
    lines += _fmt_table(generic_scores, "Top kinases — generic spectral top-1 predictions")
    lines += _fmt_table(liver_scores, "Top kinases — liver FC spectral top-1 predictions")

    # Three-way comparison
    if not psp_scores.empty and not generic_scores.empty and not liver_scores.empty:
        merged = (
            psp_scores[["kinase", "ksea_score"]]
            .rename(columns={"ksea_score": "psp_score"})
            .merge(
                generic_scores[["kinase", "ksea_score"]].rename(
                    columns={"ksea_score": "generic_score"}
                ),
                on="kinase",
                how="outer",
            )
            .merge(
                liver_scores[["kinase", "ksea_score"]].rename(
                    columns={"ksea_score": "liver_score"}
                ),
                on="kinase",
                how="outer",
            )
        )
        merged = merged.dropna().copy()
        merged["liver_vs_generic"] = merged["liver_score"] - merged["generic_score"]
        merged = merged.sort_values("psp_score", ascending=False).reset_index(drop=True)
        merged.to_csv(KSEA_DIR / "ksea_three_way_comparison.csv", index=False)

        lines += [
            "Three-way comparison (kinases scored by all three methods, sorted by PSP score):",
            "-" * 80,
            f"  {'Kinase':<20}  {'PSP score':>10}  {'Generic pred':>13}  "
            f"{'Liver pred':>11}  {'Liver-Generic':>14}",
            "  " + "-" * 74,
        ]
        for _, row in merged.head(30).iterrows():
            lines.append(
                f"  {row['kinase']:<20}  {row['psp_score']:>10.3f}  "
                f"{row['generic_score']:>13.3f}  {row['liver_score']:>11.3f}  "
                f"{row['liver_vs_generic']:>+14.3f}"
            )
        lines.append("")

        # Highlight kinases where liver FC model predicts them for more-active sites
        lines += [
            "Kinases MOST GAINED activity in liver FC model vs generic (liver_vs_generic):",
            "-" * 60,
        ]
        for _, row in merged.nlargest(10, "liver_vs_generic").iterrows():
            lines.append(
                f"  {row['kinase']:<20}  PSP={row['psp_score']:>+7.3f}  "
                f"generic={row['generic_score']:>+7.3f}  liver={row['liver_score']:>+7.3f}  "
                f"Δ={row['liver_vs_generic']:>+7.3f}"
            )
        lines.append("")

    return "\n".join(lines)


def main() -> None:
    KSEA_DIR.mkdir(parents=True, exist_ok=True)

    print("=== Loading liver FC values ===")
    fc_values = load_site_fc_values(LIVER_EXCEL)
    print(f"Phosphosites with FC values: {len(fc_values):,}")
    print(f"Global mean log2FC: {fc_values.mean():.3f}   SD: {fc_values.std():.3f}")
    print()

    print("=== Building PSP substrate sets (from generic graph) ===")
    generic_edges = read_csv_auto(GENERIC_EDGES)
    psp_sets = _build_psp_substrate_sets(generic_edges)
    print(f"Kinases with PSP substrates: {len(psp_sets):,}")

    print("=== Running KSEA — PSP substrates ===")
    psp_scores = compute_ksea_scores(fc_values, psp_sets, min_substrates=MIN_SUBSTRATES)
    psp_scores.to_csv(KSEA_DIR / "ksea_psp.csv", index=False)
    print(f"Kinases scored: {len(psp_scores):,}")
    print()

    generic_sp_path = PRED_OUTPUTS_DIR / "generic_multi_kinase_spectral" / "results.csv.gz"
    liver_sp_path = PRED_OUTPUTS_DIR / "liver_multi_kinase_spectral" / "results.csv.gz"

    generic_scores = pd.DataFrame()
    liver_scores = pd.DataFrame()

    if generic_sp_path.exists():
        print("=== Building top-1 substrate sets — generic spectral ===")
        results_g = pd.read_csv(generic_sp_path)
        pred_sets_g = _build_top1_substrate_sets(results_g)
        print(f"Kinases with top-1 predictions: {len(pred_sets_g):,}")
        print("=== Running KSEA — generic spectral top-1 predictions ===")
        generic_scores = compute_ksea_scores(fc_values, pred_sets_g, min_substrates=MIN_SUBSTRATES)
        generic_scores.to_csv(KSEA_DIR / "ksea_generic_spectral.csv", index=False)
        print(f"Kinases scored: {len(generic_scores):,}")
        print()
    else:
        print("Generic spectral results not found — skipping.")

    if liver_sp_path.exists():
        print("=== Building top-1 substrate sets — liver FC spectral ===")
        results_l = pd.read_csv(liver_sp_path)
        pred_sets_l = _build_top1_substrate_sets(results_l)
        print(f"Kinases with top-1 predictions: {len(pred_sets_l):,}")
        print("=== Running KSEA — liver FC spectral top-1 predictions ===")
        liver_scores = compute_ksea_scores(fc_values, pred_sets_l, min_substrates=MIN_SUBSTRATES)
        liver_scores.to_csv(KSEA_DIR / "ksea_liver_spectral.csv", index=False)
        print(f"Kinases scored: {len(liver_scores):,}")
        print()
    else:
        print("Liver FC spectral results not found — skipping.")

    report = _build_report(fc_values, psp_scores, generic_scores, liver_scores)
    write_text(report, KSEA_DIR / "ksea_comparison.txt")
    print(report)
    print(f"All KSEA outputs written to: {KSEA_DIR}")


if __name__ == "__main__":
    main()
