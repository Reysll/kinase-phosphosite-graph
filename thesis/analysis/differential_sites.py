"""
Differential phosphosite analysis for the cancer-center collaborator.

Approved 2026-10-25 with two changes to the original plan:
  - per-channel median normalization on log2, because
    thesis/analysis/check_normalization.py showed the Phosphorylation sheet's
    "Abundance" columns are NOT normalized as supplied (per-channel medians
    span 6.13x, 2.62 log2 units)
  - a batch-aware model as the primary test, because the design is 6 batches of
    3 tumour + 3 control and the batch effect is large and direction-varying

TESTS
-----
PRIMARY   linear model on log2:  y ~ intercept + condition + batch dummies
          OLS by numpy least squares (statsmodels is not installed and adding
          it would change the project environment). The condition coefficient
          is the batch-adjusted log2 fold change; its t statistic gives the
          p-value, with df = n - p.
SENSITIVITY  Welch unpaired t-test, ignoring batch.
SENSITIVITY  Mann-Whitney U, rank based.

The number of sites that change significance between the primary and the Welch
test is reported, since that is the direct measure of how much batch structure
mattered.

Benjamini-Hochberg FDR is applied separately within each test, over the set of
tested sites.

MISSING VALUES
--------------
No imputation. Each site is tested on the channels where it was observed.
Phosphoproteomics missingness is abundance dependent, so imputing would bias
toward the null. Sites seen in one group and never in the other cannot be
tested and get a separate presence/absence list with a Fisher exact test on
detection counts.

Inputs (read-only):
    data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx
    outputs/edges.csv.gz                      (PSP kinase-site annotation)

Writes only into thesis/analysis/out/diffsites/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.differential_sites
"""

from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

XLSX = Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")
EDGES = Path("outputs/edges.csv.gz")
OUT = Path("thesis/analysis/out/diffsites")

SITE_RE = re.compile(r"([STY])(\d+)")
CHAN_RE = re.compile(r"F(\d+)")
MIN_PER_GROUP = 6
MIN_STRICT = 12
ALPHA = 0.05
LFC_CUT = 1.0


# --------------------------------------------------------------------- utils
def bh(p: np.ndarray) -> np.ndarray:
    return stats.false_discovery_control(p, method="bh")


def hedges_g(a: np.ndarray, b: np.ndarray) -> float:
    n1, n2 = len(a), len(b)
    if n1 < 2 or n2 < 2:
        return np.nan
    s1, s2 = a.std(ddof=1), b.std(ddof=1)
    sp = np.sqrt(((n1 - 1) * s1 ** 2 + (n2 - 1) * s2 ** 2) / (n1 + n2 - 2))
    if sp == 0:
        return np.nan
    d = (b.mean() - a.mean()) / sp
    J = 1.0 - 3.0 / (4.0 * (n1 + n2) - 9.0)
    return float(d * J)


def cliffs_delta(a: np.ndarray, b: np.ndarray) -> float:
    if len(a) == 0 or len(b) == 0:
        return np.nan
    gt = (b[:, None] > a[None, :]).sum()
    lt = (b[:, None] < a[None, :]).sum()
    return float((gt - lt) / (len(a) * len(b)))


def lm_condition(y, cond, batch):
    """OLS y ~ 1 + cond + batch dummies. Returns (coef, se, t, p, df) or None."""
    ub = np.unique(batch)
    cols = [np.ones(len(y)), cond.astype(float)]
    for b in ub[1:]:
        cols.append((batch == b).astype(float))
    X = np.column_stack(cols)
    n, p = X.shape
    if n - p < 1:
        return None
    beta, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank < p:
        return None
    resid = y - X @ beta
    df = n - p
    s2 = float(resid @ resid) / df
    if s2 <= 0:
        return None
    try:
        xtx_inv = np.linalg.inv(X.T @ X)
    except np.linalg.LinAlgError:
        return None
    se = float(np.sqrt(s2 * xtx_inv[1, 1]))
    if not np.isfinite(se) or se == 0:
        return None
    t = float(beta[1]) / se
    return float(beta[1]), se, t, float(2 * stats.t.sf(abs(t), df)), int(df)


# ---------------------------------------------------------------------- main
def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    L = []
    add = L.append
    add("DIFFERENTIAL PHOSPHOSITE ANALYSIS (HCC tumour vs adjacent non-tumour)")
    add("=" * 78)
    add("")
    add("Read-only. Writes only into thesis/analysis/out/diffsites/.")
    add("")

    df = pd.read_excel(XLSX, sheet_name="Phosphorylation", header=1)
    ctrl = [c for c in df.columns if isinstance(c, str) and c.strip().endswith("Control")]
    samp = [c for c in df.columns if isinstance(c, str) and c.strip().endswith("Sample")]
    chans = ctrl + samp

    # batch id from the channel number: blocks of six consecutive F numbers
    def batch_of(col):
        m = CHAN_RE.search(col)
        return (int(m.group(1)) - 1) // 6

    batch = np.array([batch_of(c) for c in chans])
    cond = np.array([0] * len(ctrl) + [1] * len(samp))   # 0 control, 1 tumour

    add("DESIGN")
    add("-" * 78)
    add("  channels: {} control, {} tumour".format(len(ctrl), len(samp)))
    bt = pd.DataFrame({"batch": batch, "cond": cond})
    tab = bt.pivot_table(index="batch", columns="cond", aggfunc=len, fill_value=0)
    tab.columns = ["control", "tumour"]
    add("  batches inferred from channel numbering (blocks of six):")
    add("    " + tab.to_string().replace("\n", "\n    "))
    add("")

    # --------------------------------------------------- parse sites, log2, centre
    raw = df[chans].apply(pd.to_numeric, errors="coerce")
    raw = raw.where(raw > 0)
    lg = np.log2(raw)

    chan_med = lg.median(axis=0)
    lgc = lg.sub(chan_med, axis=1)

    add("NORMALIZATION")
    add("-" * 78)
    add("  Per-channel median centering applied on log2.")
    add("  Channel log2 medians before: min {:.3f}  max {:.3f}  spread {:.3f}".format(
        chan_med.min(), chan_med.max(), chan_med.max() - chan_med.min()))
    after = lgc.median(axis=0)
    add("  Channel log2 medians after : min {:.3f}  max {:.3f}  spread {:.3f}".format(
        after.min(), after.max(), after.max() - after.min()))
    add("  See normalization_check.txt for the diagnostic that motivated this.")
    add("")

    nsite = df["Modifications in Master Proteins"].apply(
        lambda v: len(SITE_RE.findall(v)) if isinstance(v, str) else 0)
    one = nsite == 1
    gene = df["Gene Symbol"].astype(str).str.strip()
    lab = df["Modifications in Master Proteins"].apply(
        lambda v: "".join(SITE_RE.findall(v)[0]) if isinstance(v, str)
        and len(SITE_RE.findall(v)) == 1 else None)
    site_id = gene + "-" + lab.astype(str)

    add("SITE PARSING")
    add("-" * 78)
    add("  rows in sheet                         : {:,}".format(len(df)))
    add("  rows with exactly one phosphosite     : {:,}".format(int(one.sum())))
    add("  rows with more than one (excluded)    : {:,}".format(int((nsite > 1).sum())))
    add("  rows with none parseable (excluded)   : {:,}".format(int((nsite == 0).sum())))

    sub = lgc.loc[one].copy()
    sub["site"] = site_id[one].to_numpy()
    collapsed = sub.groupby("site", observed=True)[chans].mean()
    add("  unique sites after collapsing peptides: {:,}".format(len(collapsed)))
    add("    duplicate peptide rows are averaged per channel on the log2 scale,")
    add("    so each site gets exactly one test and the FDR denominator is sites")
    add("")

    M = collapsed[chans].to_numpy()
    is_c = cond == 0
    is_t = cond == 1
    n_c = np.isfinite(M[:, is_c]).sum(axis=1)
    n_t = np.isfinite(M[:, is_t]).sum(axis=1)

    testable = (n_c >= MIN_PER_GROUP) & (n_t >= MIN_PER_GROUP)
    add("TESTABLE SET")
    add("-" * 78)
    add("  sites with >= {} observed in BOTH groups: {:,}".format(
        MIN_PER_GROUP, int(testable.sum())))
    add("  sites with >= {} observed in BOTH groups: {:,}  (strict tier)".format(
        MIN_STRICT, int(((n_c >= MIN_STRICT) & (n_t >= MIN_STRICT)).sum())))
    add("")

    # ------------------------------------------------------------------ tests
    idx = np.where(testable)[0]
    names = collapsed.index.to_numpy()
    recs = []
    for i in idx:
        row = M[i]
        obs = np.isfinite(row)
        y = row[obs]
        c = cond[obs]
        b = batch[obs]
        a_ = y[c == 0]
        b_ = y[c == 1]

        lm = lm_condition(y, c, b)
        tw, pw = stats.ttest_ind(b_, a_, equal_var=False)
        try:
            mu, pm = stats.mannwhitneyu(b_, a_, alternative="two-sided")
        except ValueError:
            mu, pm = np.nan, np.nan

        recs.append({
            "site": names[i],
            "gene": str(names[i]).rsplit("-", 1)[0],
            "residue": str(names[i]).rsplit("-", 1)[-1],
            "n_control": int(a_.size),
            "n_tumour": int(b_.size),
            "n_batches": int(len(np.unique(b))),
            "mean_control_log2": float(a_.mean()),
            "mean_tumour_log2": float(b_.mean()),
            "log2fc_raw": float(b_.mean() - a_.mean()),
            "log2fc_adj": lm[0] if lm else np.nan,
            "lm_se": lm[1] if lm else np.nan,
            "lm_t": lm[2] if lm else np.nan,
            "lm_p": lm[3] if lm else np.nan,
            "lm_df": lm[4] if lm else np.nan,
            "welch_t": float(tw),
            "welch_p": float(pw),
            "mwu_p": float(pm) if pm == pm else np.nan,
            "hedges_g": hedges_g(a_, b_),
            "cliffs_delta": cliffs_delta(a_, b_),
        })

    res = pd.DataFrame(recs)
    for src, dst in (("lm_p", "lm_fdr"), ("welch_p", "welch_fdr"), ("mwu_p", "mwu_fdr")):
        ok = res[src].notna().to_numpy()
        res[dst] = np.nan
        res.loc[ok, dst] = bh(res.loc[ok, src].to_numpy())

    res["sig_lm"] = res["lm_fdr"] < ALPHA
    res["sig_welch"] = res["welch_fdr"] < ALPHA
    res["sig_mwu"] = res["mwu_fdr"] < ALPHA
    res["shortlist"] = (res["sig_lm"] & (res["log2fc_adj"].abs() >= LFC_CUT)
                        & (res["n_control"] >= 9) & (res["n_tumour"] >= 9))

    add("RESULTS: PRIMARY TEST (linear model, condition + batch)")
    add("-" * 78)
    add("  sites tested                 : {:,}".format(len(res)))
    add("  model not estimable          : {:,}".format(int(res["lm_p"].isna().sum())))
    add("  significant at BH-FDR < 0.05 : {:,}  ({:.1f}%)".format(
        int(res["sig_lm"].sum()), 100 * res["sig_lm"].mean()))
    up = int((res["sig_lm"] & (res["log2fc_adj"] > 0)).sum())
    dn = int((res["sig_lm"] & (res["log2fc_adj"] < 0)).sum())
    add("    up in tumour {:,}   down in tumour {:,}".format(up, dn))
    add("  |log2FC| >= 1 among significant : {:,}".format(
        int((res["sig_lm"] & (res["log2fc_adj"].abs() >= 1)).sum())))
    add("  shortlist (FDR<0.05, |log2FC|>=1, n>=9 both groups): {:,}".format(
        int(res["shortlist"].sum())))
    add("")

    add("RESULTS: SENSITIVITY TESTS")
    add("-" * 78)
    add("  Welch unpaired, ignoring batch : {:,} significant ({:.1f}%)".format(
        int(res["sig_welch"].sum()), 100 * res["sig_welch"].mean()))
    add("  Mann-Whitney U                 : {:,} significant ({:.1f}%)".format(
        int(res["sig_mwu"].sum()), 100 * res["sig_mwu"].mean()))
    add("")

    both = int((res["sig_lm"] & res["sig_welch"]).sum())
    lm_only = int((res["sig_lm"] & ~res["sig_welch"]).sum())
    w_only = int((~res["sig_lm"] & res["sig_welch"]).sum())
    neither = int((~res["sig_lm"] & ~res["sig_welch"]).sum())
    add("HOW MUCH DID BATCH STRUCTURE MATTER?")
    add("-" * 78)
    add("  significant in BOTH                  : {:,}".format(both))
    add("  gained by modelling batch (LM only)  : {:,}".format(lm_only))
    add("  lost by modelling batch (Welch only) : {:,}".format(w_only))
    add("  significant in neither                : {:,}".format(neither))
    add("  sites that CHANGE significance        : {:,}  ({:.1f}% of tested)".format(
        lm_only + w_only, 100.0 * (lm_only + w_only) / len(res)))
    add("")
    add("  correlation of log2FC, raw vs batch-adjusted: r = {:.4f}".format(
        res["log2fc_raw"].corr(res["log2fc_adj"])))
    add("")

    sw = res[res["sig_lm"] & ~res["sig_welch"]].nsmallest(10, "lm_fdr")
    if not sw.empty:
        add("  10 sites most clearly rescued by the batch model:")
        d = sw[["site", "n_control", "n_tumour", "log2fc_raw", "log2fc_adj",
                "lm_fdr", "welch_fdr"]].copy()
        for c in ("log2fc_raw", "log2fc_adj"):
            d[c] = d[c].map("{:+.3f}".format)
        for c in ("lm_fdr", "welch_fdr"):
            d[c] = d[c].map(lambda v: "{:.2e}".format(v))
        add("    " + d.to_string(index=False).replace("\n", "\n    "))
        add("")

    add("TOP 30 SITES BY PRIMARY TEST")
    add("-" * 78)
    top = res.nsmallest(30, "lm_fdr")[
        ["site", "n_control", "n_tumour", "log2fc_adj", "hedges_g",
         "cliffs_delta", "lm_fdr", "sig_welch"]].copy()
    top["log2fc_adj"] = top["log2fc_adj"].map("{:+.3f}".format)
    top["hedges_g"] = top["hedges_g"].map("{:+.2f}".format)
    top["cliffs_delta"] = top["cliffs_delta"].map("{:+.2f}".format)
    top["lm_fdr"] = top["lm_fdr"].map(lambda v: "{:.2e}".format(v))
    add(top.to_string(index=False))
    add("")

    # ------------------------------------------------- presence / absence
    pa_mask = ((n_c == 0) & (n_t >= MIN_PER_GROUP)) | ((n_t == 0) & (n_c >= MIN_PER_GROUP))
    pa_rows = []
    for i in np.where(pa_mask)[0]:
        table = [[int(n_t[i]), int(len(samp) - n_t[i])],
                 [int(n_c[i]), int(len(ctrl) - n_c[i])]]
        _, pv = stats.fisher_exact(table)
        pa_rows.append({
            "site": names[i], "gene": str(names[i]).rsplit("-", 1)[0],
            "n_control": int(n_c[i]), "n_tumour": int(n_t[i]),
            "direction": "tumour only" if n_c[i] == 0 else "control only",
            "fisher_p": float(pv),
        })
    pa = pd.DataFrame(pa_rows)
    if not pa.empty:
        pa["fisher_fdr"] = bh(pa["fisher_p"].to_numpy())
        pa["significant"] = pa["fisher_fdr"] < ALPHA
        pa = pa.sort_values("fisher_fdr")
    add("PRESENCE / ABSENCE (detected in one group only, never tested above)")
    add("-" * 78)
    add("  sites: {:,}".format(len(pa)))
    if not pa.empty:
        add("    tumour only  : {:,}".format(int((pa["direction"] == "tumour only").sum())))
        add("    control only : {:,}".format(int((pa["direction"] == "control only").sum())))
        add("    significant at BH-FDR < 0.05: {:,}".format(int(pa["significant"].sum())))
        add("")
        d = pa.head(15)[["site", "n_control", "n_tumour", "direction", "fisher_fdr"]].copy()
        d["fisher_fdr"] = d["fisher_fdr"].map(lambda v: "{:.2e}".format(v))
        add("  " + d.to_string(index=False).replace("\n", "\n  "))
    add("")

    # --------------------------------------------------- kinase mapping (PSP)
    edges = pd.read_csv(EDGES, usecols=["source", "target", "relation"])
    psp = edges[edges["relation"] == "phosphorylates"]
    site2kin = psp.groupby("target")["source"].apply(
        lambda s: sorted({x.split(":")[-1] for x in s})).to_dict()

    res["psp_site_id"] = "SITE:" + res["site"]
    res["psp_kinases"] = res["psp_site_id"].map(
        lambda s: ";".join(site2kin.get(s, [])))
    res["n_psp_kinases"] = res["psp_site_id"].map(lambda s: len(site2kin.get(s, [])))

    sig = res[res["sig_lm"]]
    add("KINASE MAPPING, ROUTE A: DIRECT PhosphoSitePlus ANNOTATION")
    add("-" * 78)
    add("  tested sites with >=1 annotated kinase      : {:,} of {:,}  ({:.1f}%)".format(
        int((res["n_psp_kinases"] > 0).sum()), len(res),
        100 * (res["n_psp_kinases"] > 0).mean()))
    add("  SIGNIFICANT sites with >=1 annotated kinase : {:,} of {:,}  ({:.1f}%)".format(
        int((sig["n_psp_kinases"] > 0).sum()), len(sig),
        100 * (sig["n_psp_kinases"] > 0).mean() if len(sig) else 0.0))
    add("  Coverage is stated rather than filtered: unannotated sites remain in")
    add("  every output, since the collaborator may still want to assay them.")
    kin_counts = {}
    for s in sig["psp_kinases"]:
        for k in [x for x in s.split(";") if x]:
            kin_counts[k] = kin_counts.get(k, 0) + 1
    if kin_counts:
        kc = pd.Series(kin_counts).sort_values(ascending=False)
        add("")
        add("  kinases with the most significant annotated sites:")
        add("    " + kc.head(20).to_string().replace("\n", "\n    "))
    add("")

    # ---------------------------------------------------- KSEA (two-tailed)
    lfc = res.set_index("psp_site_id")["log2fc_adj"]
    gmean = float(lfc.mean())
    delta = float(lfc.std(ddof=1))
    krows = []
    kin2sites = {}
    for site, kins in site2kin.items():
        if site in lfc.index:
            for k in kins:
                kin2sites.setdefault(k, []).append(site)
    for k, sites in kin2sites.items():
        vals = lfc.loc[sites].dropna()
        m = len(vals)
        if m < 3:
            continue
        score = (vals.mean() - gmean) * np.sqrt(m) / delta
        krows.append({
            "kinase": k, "n_substrates_measured": m,
            "mean_substrate_log2fc": float(vals.mean()),
            "global_mean_log2fc": gmean,
            "ksea_score": float(score),
            "p_value": float(2 * stats.norm.sf(abs(score))),
        })
    ksea = pd.DataFrame(krows)
    if not ksea.empty:
        ksea["adj_p"] = bh(ksea["p_value"].to_numpy())
        ksea["significant"] = ksea["adj_p"] < ALPHA
        ksea = ksea.sort_values("ksea_score", ascending=False)
    add("KINASE MAPPING, ROUTE B: KSEA ON THE DIFFERENTIAL RESULT (two-tailed)")
    add("-" * 78)
    add("  score = (mean substrate log2FC - global mean) * sqrt(m) / sd")
    add("  p = 2 * norm.sf(|score|). The published 17/64 used a one-tailed p.")
    add("  kinases scoreable (>=3 measured substrates): {:,}".format(len(ksea)))
    if not ksea.empty:
        add("  significant at BH-FDR < 0.05             : {:,}".format(
            int(ksea["significant"].sum())))
        add("")
        d = ksea.head(12)[["kinase", "n_substrates_measured",
                           "mean_substrate_log2fc", "ksea_score", "adj_p"]].copy()
        d["mean_substrate_log2fc"] = d["mean_substrate_log2fc"].map("{:+.3f}".format)
        d["ksea_score"] = d["ksea_score"].map("{:+.2f}".format)
        d["adj_p"] = d["adj_p"].map(lambda v: "{:.2e}".format(v))
        add("  most UP in tumour:")
        add("    " + d.to_string(index=False).replace("\n", "\n    "))
        d = ksea.tail(12).iloc[::-1][["kinase", "n_substrates_measured",
                                      "mean_substrate_log2fc", "ksea_score", "adj_p"]].copy()
        d["mean_substrate_log2fc"] = d["mean_substrate_log2fc"].map("{:+.3f}".format)
        d["ksea_score"] = d["ksea_score"].map("{:+.2f}".format)
        d["adj_p"] = d["adj_p"].map(lambda v: "{:.2e}".format(v))
        add("")
        add("  most DOWN in tumour:")
        add("    " + d.to_string(index=False).replace("\n", "\n    "))
    add("")

    # ------------------------------- disjointness of the two kinase routes
    add("DO THE TWO ROUTES NOMINATE THE SAME KINASES?")
    add("-" * 78)
    route_a = set(kin_counts.keys())
    route_b = set(ksea.loc[ksea["significant"], "kinase"]) if not ksea.empty else set()
    inter = route_a & route_b
    union = route_a | route_b
    add("  Route A, kinases annotated to >=1 significant site : {:,}".format(len(route_a)))
    add("  Route B, KSEA-significant kinases                  : {:,}".format(len(route_b)))
    add("  overlap                                            : {:,}".format(len(inter)))
    if union:
        add("  Jaccard                                            : {:.3f}".format(
            len(inter) / len(union)))
    if route_b:
        add("  share of Route B found in Route A                  : {:.1f}%".format(
            100.0 * len(inter) / len(route_b)))
    if inter:
        add("  overlapping kinases: " + ", ".join(sorted(inter)[:30]))
    add("")

    # ------------------------------------------------------------- volcano
    try:
        fig, ax = plt.subplots(figsize=(9, 7), dpi=150)
        x = res["log2fc_adj"].to_numpy()
        y = -np.log10(res["lm_fdr"].to_numpy())
        ns = ~res["sig_lm"].to_numpy()
        ax.scatter(x[ns], y[ns], s=8, c="#BBBBBB", linewidths=0, label="not significant")
        s_ = res["sig_lm"].to_numpy() & ~res["shortlist"].to_numpy()
        ax.scatter(x[s_], y[s_], s=10, c="#4C72B0", linewidths=0, label="FDR < 0.05")
        sl = res["shortlist"].to_numpy()
        ax.scatter(x[sl], y[sl], s=22, c="#C24E10", linewidths=0,
                   label="shortlist (|log2FC| >= 1, n >= 9)")
        ax.axhline(-np.log10(ALPHA), ls="--", lw=0.8, c="#777777")
        ax.axvline(LFC_CUT, ls="--", lw=0.8, c="#777777")
        ax.axvline(-LFC_CUT, ls="--", lw=0.8, c="#777777")
        for _, r in res[res["shortlist"]].nsmallest(15, "lm_fdr").iterrows():
            ax.annotate(r["site"], (r["log2fc_adj"], -np.log10(r["lm_fdr"])),
                        fontsize=6, alpha=0.85,
                        xytext=(3, 3), textcoords="offset points")
        ax.set_xlabel("batch-adjusted log2 fold change (tumour / adjacent normal)")
        ax.set_ylabel("-log10 BH-FDR")
        ax.set_title("Differential phosphosites, HCC tumour vs adjacent non-tumour\n"
                     "linear model with batch, {:,} sites tested".format(len(res)))
        ax.legend(fontsize=8, loc="lower left", framealpha=0.9)
        fig.tight_layout()
        fig.savefig(OUT / "diffsites_volcano.png")
        plt.close(fig)
        add("Volcano written: diffsites_volcano.png")
    except Exception as e:                                     # noqa: BLE001
        add("Volcano FAILED: {}".format(e))
    add("")

    # --------------------------------------------------------------- outputs
    strict = (n_c >= MIN_STRICT) & (n_t >= MIN_STRICT)
    res["strict_tier"] = res["site"].isin(names[strict])

    res.to_csv(OUT / "diffsites_tests.csv", index=False)
    res[res["sig_lm"]].sort_values("lm_fdr").to_csv(
        OUT / "diffsites_significant.csv", index=False)
    res[["site", "gene", "residue", "log2fc_adj", "lm_fdr", "sig_lm", "sig_welch",
         "sig_mwu", "strict_tier"]].to_csv(OUT / "diffsites_sensitivity.csv", index=False)
    if not pa.empty:
        pa.to_csv(OUT / "diffsites_presence_absence.csv", index=False)
    res.loc[res["n_psp_kinases"] > 0,
            ["site", "gene", "residue", "log2fc_adj", "lm_fdr", "sig_lm",
             "psp_kinases", "n_psp_kinases"]].to_csv(
        OUT / "diffsites_kinase_map.csv", index=False)
    if not ksea.empty:
        ksea.to_csv(OUT / "diffsites_ksea.csv", index=False)

    add("FILES WRITTEN (all under thesis/analysis/out/diffsites/)")
    add("-" * 78)
    for f in ("diffsites.txt", "diffsites_tests.csv", "diffsites_significant.csv",
              "diffsites_sensitivity.csv", "diffsites_presence_absence.csv",
              "diffsites_kinase_map.csv", "diffsites_ksea.csv",
              "diffsites_volcano.png", "normalization_check.txt",
              "channel_stats.csv"):
        add("  " + f)

    text = "\n".join(L)
    (OUT / "diffsites.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
