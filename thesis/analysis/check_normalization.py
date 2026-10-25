"""
Step 3 of the differential-site plan: is the Phosphorylation sheet normalized?

The ProteinExpression sheet's columns are labelled "Abundances (Normalized)".
The Phosphorylation sheet's are labelled only "Abundance". That difference in
labelling may or may not reflect a difference in processing, so this script
measures it rather than assuming.

Diagnostics, per sheet:
  - per-channel median, mean, IQR and total signal on the raw scale
  - the same on log2, which is the scale any test would use
  - spread of per-channel medians across channels (max/min ratio, CV)
  - per-channel count of observed values

If a sheet is median-normalized, its per-channel medians should be nearly
identical. If they are not, per-channel median centering on log2 is needed
before any across-channel comparison.

Read-only. Writes only into thesis/analysis/out/diffsites/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.check_normalization
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

XLSX = Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")
OUT = Path("thesis/analysis/out/diffsites")

SHEETS = {
    "ProteinExpression": "Abundances (Normalized)",
    "Phosphorylation": "Abundance",
}


def split_cols(cols):
    ctrl = [c for c in cols if isinstance(c, str) and c.strip().endswith("Control")]
    samp = [c for c in cols if isinstance(c, str) and c.strip().endswith("Sample")]
    return ctrl, samp


def channel_stats(mat: pd.DataFrame) -> pd.DataFrame:
    rows = []
    lg = np.log2(mat.where(mat > 0))
    for c in mat.columns:
        v = mat[c].dropna()
        l = lg[c].dropna()
        rows.append({
            "channel": c,
            "n_observed": int(v.size),
            "median": float(v.median()) if v.size else np.nan,
            "mean": float(v.mean()) if v.size else np.nan,
            "q1": float(v.quantile(0.25)) if v.size else np.nan,
            "q3": float(v.quantile(0.75)) if v.size else np.nan,
            "total_signal": float(v.sum()) if v.size else np.nan,
            "log2_median": float(l.median()) if l.size else np.nan,
            "log2_sd": float(l.std()) if l.size else np.nan,
        })
    return pd.DataFrame(rows)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L = []
    add = L.append
    add("ARE THE PHOSPHORYLATION ABUNDANCES NORMALIZED?")
    add("=" * 78)
    add("")
    add("Read-only diagnostic. No test is run here.")
    add("")

    verdicts = {}
    frames = {}

    for sheet, label in SHEETS.items():
        df = pd.read_excel(XLSX, sheet_name=sheet, header=1)
        ctrl, samp = split_cols(df.columns)
        cols = ctrl + samp
        mat = df[cols].apply(pd.to_numeric, errors="coerce")
        st = channel_stats(mat)
        st["group"] = ["control"] * len(ctrl) + ["tumour"] * len(samp)
        st["sheet"] = sheet
        frames[sheet] = st

        med = st["median"].to_numpy()
        lmed = st["log2_median"].to_numpy()
        ratio = float(np.nanmax(med) / np.nanmin(med))
        cv = float(np.nanstd(med) / np.nanmean(med))
        lspread = float(np.nanmax(lmed) - np.nanmin(lmed))

        add("-" * 78)
        add("SHEET: {}   column label: \"{}\"".format(sheet, label))
        add("-" * 78)
        add("  channels: {} control + {} tumour".format(len(ctrl), len(samp)))
        add("  rows: {:,}".format(len(df)))
        add("  observed values per channel: min {:,}  median {:,}  max {:,}".format(
            int(st["n_observed"].min()), int(st["n_observed"].median()),
            int(st["n_observed"].max())))
        add("")
        add("  per-channel MEDIAN on the raw scale:")
        add("    min {:.4g}   max {:.4g}   max/min = {:.3f}".format(
            np.nanmin(med), np.nanmax(med), ratio))
        add("    coefficient of variation across channels = {:.4f}".format(cv))
        add("  per-channel MEDIAN on log2:")
        add("    spread (max - min) = {:.4f} log2 units".format(lspread))
        add("")

        normalized = ratio < 1.10
        verdicts[sheet] = {
            "ratio": ratio, "cv": cv, "log2_spread": lspread,
            "looks_normalized": normalized,
        }
        add("  VERDICT: per-channel medians {} (max/min = {:.3f}).".format(
            "agree to within 10%" if normalized else "DIFFER by more than 10%",
            ratio))
        add("           -> {}".format(
            "consistent with median normalization already applied"
            if normalized else
            "NOT median-normalized; per-channel centering is required"))
        add("")
        add("  per-channel detail:")
        disp = st[["channel", "group", "n_observed", "median", "q1", "q3",
                   "log2_median", "log2_sd"]].copy()
        disp["channel"] = disp["channel"].str.replace("Abundances (Normalized): ", "",
                                                      regex=False)
        disp["channel"] = disp["channel"].str.replace("Abundance: ", "", regex=False)
        for c in ("median", "q1", "q3"):
            disp[c] = disp[c].map(lambda v: "{:.4g}".format(v))
        for c in ("log2_median", "log2_sd"):
            disp[c] = disp[c].map(lambda v: "{:.3f}".format(v))
        add(disp.to_string(index=False))
        add("")

    allst = pd.concat(frames.values(), ignore_index=True)
    allst.to_csv(OUT / "channel_stats.csv", index=False)

    add("=" * 78)
    add("CONCLUSION")
    add("=" * 78)
    ph = verdicts["Phosphorylation"]
    pr = verdicts["ProteinExpression"]
    add("ProteinExpression, labelled normalized : max/min median = {:.3f}".format(
        pr["ratio"]))
    add("Phosphorylation,  labelled abundance   : max/min median = {:.3f}".format(
        ph["ratio"]))
    add("")
    if ph["looks_normalized"]:
        add("The phosphorylation channels are already comparable. Per-channel")
        add("median centering would be a no-op, but it is harmless and will be")
        add("applied anyway so the test does not depend on this judgement.")
    else:
        add("The phosphorylation channels are NOT comparable as supplied.")
        add("Per-channel median centering on log2 WILL be applied before testing:")
        add("    log2_centered[:, j] = log2(x[:, j]) - median(log2(x[:, j]))")
        add("computed over the sites observed in channel j. This removes")
        add("channel loading differences without touching within-channel ranking.")
    add("")
    add("Written: thesis/analysis/out/diffsites/channel_stats.csv")

    text = "\n".join(L)
    (OUT / "normalization_check.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
