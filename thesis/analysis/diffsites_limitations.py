"""
Limitation checks for the differential phosphosite analysis (read-only).

Two quantities requested after review:

1. Protein-abundance confounding. A phosphosite can change because the
   underlying protein changed, not because its phosphorylation stoichiometry
   changed. Correcting for this requires protein-level abundance for the same
   samples. The ProteinExpression sheet cannot serve: its caption describes it
   as proteins unique to tumour tissue, it holds 444 rows, and most cells are
   empty. This measures exactly how many shortlist sites have their protein
   present in that sheet at all.

2. Presence/absence batch structure. A site called "tumour only" is more
   convincing if the control channels in those same batches were producing
   measurements for other sites at the time, and less convincing if the
   tumour detections sit in batches where control coverage was poor. This
   reports, for each presence/absence site, how its detections distribute
   across the six batches and what the opposite group's coverage was there.

Inputs (read-only):
    data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx
    thesis/analysis/out/diffsites/diffsites_tests.csv
    thesis/analysis/out/diffsites/diffsites_presence_absence.csv

Writes only into thesis/analysis/out/diffsites/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.diffsites_limitations
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

XLSX = Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")
OUT = Path("thesis/analysis/out/diffsites")
SITE_RE = re.compile(r"([STY])(\d+)")
CHAN_RE = re.compile(r"F(\d+)")


def main() -> None:
    L = []
    add = L.append
    add("DIFFERENTIAL ANALYSIS: LIMITATION CHECKS")
    add("=" * 78)
    add("")

    res = pd.read_csv(OUT / "diffsites_tests.csv")
    short = res[res["shortlist"]]

    # ------------------------------------------------ 1. protein coverage
    pe = pd.read_excel(XLSX, sheet_name="ProteinExpression", header=1)
    pe_ctrl = [c for c in pe.columns if isinstance(c, str) and c.strip().endswith("Control")]
    pe_samp = [c for c in pe.columns if isinstance(c, str) and c.strip().endswith("Sample")]
    pe_mat = pe[pe_ctrl + pe_samp].apply(pd.to_numeric, errors="coerce")
    pe_genes = set(pe["Gene Symbol"].astype(str).str.strip())
    pe_obs = pe_mat.notna().sum(axis=1)
    pe["n_obs"] = pe_obs

    add("1. PROTEIN-ABUNDANCE CONFOUNDING")
    add("-" * 78)
    add("  The ProteinExpression sheet:")
    add("    rows                                   : {:,}".format(len(pe)))
    add("    distinct gene symbols                  : {:,}".format(len(pe_genes)))
    add("    channels                               : {} control + {} tumour".format(
        len(pe_ctrl), len(pe_samp)))
    add("    observed values per row: median {:.0f} of {}, mean {:.1f}".format(
        pe_obs.median(), len(pe_ctrl) + len(pe_samp), pe_obs.mean()))
    add("    rows with >= 6 observed in BOTH groups : {:,}".format(int((
        (pe_mat[pe_ctrl].notna().sum(axis=1) >= 6) &
        (pe_mat[pe_samp].notna().sum(axis=1) >= 6)).sum())))
    add("    overall missingness                    : {:.1f}%".format(
        100 * (1 - pe_mat.notna().to_numpy().sum() / pe_mat.size)))
    add("")

    for label, frame in (("all tested sites", res), ("shortlist sites", short)):
        genes = frame["gene"].astype(str)
        present = genes.isin(pe_genes)
        add("  {}: {:,}".format(label, len(frame)))
        add("    distinct genes                       : {:,}".format(genes.nunique()))
        add("    sites whose protein appears in the sheet at all : {:,}  ({:.1f}%)".format(
            int(present.sum()), 100 * present.mean()))
        add("    distinct genes present               : {:,}".format(
            genes[present].nunique()))
        add("")

    # how many shortlist genes are usable for a correction
    usable = pe.loc[
        (pe_mat[pe_ctrl].notna().sum(axis=1) >= 6) &
        (pe_mat[pe_samp].notna().sum(axis=1) >= 6), "Gene Symbol"
    ].astype(str).str.strip()
    usable = set(usable)
    sg = short["gene"].astype(str)
    add("  shortlist sites whose protein is quantified well enough to support")
    add("  a correction (>= 6 observed in both groups): {:,} of {:,}  ({:.1f}%)".format(
        int(sg.isin(usable).sum()), len(short), 100 * sg.isin(usable).mean()))
    add("")
    add("  CONCLUSION: phosphosite changes reported here are NOT corrected for")
    add("  protein abundance. A site may change because its protein changed.")
    add("  The ProteinExpression sheet cannot support a correction, both because")
    add("  of coverage and because its caption describes it as proteins unique to")
    add("  tumour tissue rather than a matched quantitative matrix. Any site taken")
    add("  forward for validation should have its protein level checked")
    add("  independently.")
    add("")

    # ------------------------------------- 2. presence/absence batch structure
    pa_path = OUT / "diffsites_presence_absence.csv"
    add("2. PRESENCE/ABSENCE SITES: BATCH STRUCTURE")
    add("-" * 78)
    if not pa_path.exists():
        add("  no presence/absence file found")
    else:
        pa = pd.read_csv(pa_path)
        ph = pd.read_excel(XLSX, sheet_name="Phosphorylation", header=1)
        ctrl = [c for c in ph.columns if isinstance(c, str) and c.strip().endswith("Control")]
        samp = [c for c in ph.columns if isinstance(c, str) and c.strip().endswith("Sample")]
        chans = ctrl + samp
        batch = {c: (int(CHAN_RE.search(c).group(1)) - 1) // 6 for c in chans}
        batches = sorted(set(batch.values()))

        mat = ph[chans].apply(pd.to_numeric, errors="coerce")
        mat = mat.where(mat > 0)
        nsite = ph["Modifications in Master Proteins"].apply(
            lambda v: len(SITE_RE.findall(v)) if isinstance(v, str) else 0)
        one = nsite == 1
        lab = ph["Modifications in Master Proteins"].apply(
            lambda v: "".join(SITE_RE.findall(v)[0]) if isinstance(v, str)
            and len(SITE_RE.findall(v)) == 1 else None)
        sid = ph["Gene Symbol"].astype(str).str.strip() + "-" + lab.astype(str)
        sub = mat.loc[one].copy()
        sub["site"] = sid[one].to_numpy()
        coll = sub.groupby("site", observed=True)[chans].mean()

        # per-batch coverage of the whole dataset, for reference
        add("  Dataset-wide coverage per batch (share of all {:,} collapsed sites".format(
            len(coll)))
        add("  observed in at least one channel of that batch and group):")
        ref = []
        for b in batches:
            bc = [c for c in ctrl if batch[c] == b]
            bs = [c for c in samp if batch[c] == b]
            ref.append({
                "batch": b,
                "control_cov_pct": 100 * coll[bc].notna().any(axis=1).mean(),
                "tumour_cov_pct": 100 * coll[bs].notna().any(axis=1).mean(),
            })
        refdf = pd.DataFrame(ref)
        add("    " + refdf.round(1).to_string(index=False).replace("\n", "\n    "))
        add("")

        rows = []
        for _, r in pa.iterrows():
            site = r["site"]
            if site not in coll.index:
                continue
            v = coll.loc[site]
            det_group = samp if r["direction"] == "tumour only" else ctrl
            opp_group = ctrl if r["direction"] == "tumour only" else samp
            nb = 0
            opp_ok = 0
            for b in batches:
                dchan = [c for c in det_group if batch[c] == b]
                if v[dchan].notna().any():
                    nb += 1
                    ochan = [c for c in opp_group if batch[c] == b]
                    # was the opposite group producing data in that batch at all?
                    if coll[ochan].notna().any(axis=1).mean() > 0.5:
                        opp_ok += 1
            rows.append({
                "site": site, "direction": r["direction"],
                "n_control": r["n_control"], "n_tumour": r["n_tumour"],
                "batches_with_detection": nb,
                "of_those_opposite_group_well_covered": opp_ok,
            })
        bt = pd.DataFrame(rows)
        if not bt.empty:
            bt.to_csv(OUT / "diffsites_presence_absence_batches.csv", index=False)
            add("  presence/absence sites examined : {:,}".format(len(bt)))
            add("  detections spanning all 6 batches: {:,}  ({:.1f}%)".format(
                int((bt["batches_with_detection"] == 6).sum()),
                100 * (bt["batches_with_detection"] == 6).mean()))
            add("  detections spanning >= 4 batches : {:,}  ({:.1f}%)".format(
                int((bt["batches_with_detection"] >= 4).sum()),
                100 * (bt["batches_with_detection"] >= 4).mean()))
            add("  detections confined to 1 batch   : {:,}  ({:.1f}%)".format(
                int((bt["batches_with_detection"] == 1).sum()),
                100 * (bt["batches_with_detection"] == 1).mean()))
            add("")
            agree = (bt["batches_with_detection"] ==
                     bt["of_those_opposite_group_well_covered"])
            add("  sites where EVERY batch containing a detection also had the")
            add("  opposite group well covered: {:,}  ({:.1f}%)".format(
                int(agree.sum()), 100 * agree.mean()))
            add("")
            add("  distribution of batches with detection:")
            vc = bt["batches_with_detection"].value_counts().sort_index()
            for k, v in vc.items():
                add("    {} batch(es): {:,} sites".format(k, int(v)))
            add("")
        add("  CONCLUSION: these 114 sites are selected BY a detection-count")
        add("  filter and then tested WITH a detection-count test, so their")
        add("  significance is largely a property of the filter. They are reported")
        add("  as candidates for follow-up, not as a statistical finding.")

    text = "\n".join(L)
    (OUT / "diffsites_limitations.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
