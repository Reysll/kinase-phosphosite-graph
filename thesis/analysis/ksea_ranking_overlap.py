"""
Do KSEA and the ranking model nominate the same kinases? (read-only)

STATE.md open decision 6 flagged a sentence in Chapter II section 2.7 asserting
that "KSEA and ranking nominate largely disjoint kinase sets". The assertion was
plausible from the project's KSEA report and the top-1 tables but had never been
computed. This script computes it.

Two nominated sets are compared, both restricted to the kinases each method
could in principle nominate:

  KSEA set      kinases significant at BH-FDR < 0.05 using the TWO-TAILED
                p-value, 2 * norm.sf(|score|). The project's published 17 of 64
                used a one-tailed p; the two-tailed count is 12.

  Ranking set   kinases the predict-all pass ranks highly. Reported at three
                cutoffs, since "nominate" has no single definition: the kinases
                whose median rank_in_site falls in the top 10, top 20 and top 50
                of the 415 scored.

The comparison universe is the intersection of the two methods' scoreable
kinases, because a kinase KSEA cannot score (too few measured substrates) is not
evidence of disagreement. Overlap is tested against chance with a hypergeometric
test, two-sided plus the lower tail, the lower tail being the one that would
support a claim of disjointness.

Inputs (read-only):
    thesis/analysis/out/ksea_twotailed.csv
    outputs_prediction/inference_all/{control,liver_fc}_spectral/ranks.csv.gz

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.ksea_ranking_overlap
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUT = Path("thesis/analysis/out")
INFER = Path("outputs_prediction/inference_all")
KSEA = OUT / "ksea_twotailed.csv"
CUTOFFS = (10, 20, 50)


def main() -> None:
    L = []
    add = L.append
    add("KSEA vs RANKING: DO THEY NOMINATE THE SAME KINASES?")
    add("=" * 78)
    add("")

    ks = pd.read_csv(KSEA)
    ks["kinase"] = ks["kinase"].astype(str)
    ksea_scoreable = set(ks["kinase"])
    ksea_sig = set(ks.loc[ks["adj_p_two_tailed"] < 0.05, "kinase"])
    ksea_sig_1t = set(ks.loc[ks["adj_p_recomputed_one_tailed"] < 0.05, "kinase"])

    add("KSEA (two-tailed, BH-FDR < 0.05)")
    add("-" * 78)
    add("  scoreable kinases          : {}".format(len(ksea_scoreable)))
    add("  significant, two-tailed    : {}".format(len(ksea_sig)))
    add("  significant, one-tailed    : {}  (the figure in older project docs)".format(
        len(ksea_sig_1t)))
    add("  two-tailed set: " + ", ".join(sorted(ksea_sig)))
    add("")

    rows = []
    for net in ("control", "liver_fc"):
        path = INFER / "{}_spectral".format(net) / "ranks.csv.gz"
        if not path.exists():
            add("MISSING: {}".format(path))
            continue
        df = pd.read_csv(path, usecols=["kinase_node_id", "rank_in_site"])
        df["kinase"] = df["kinase_node_id"].astype(str).str.split(":").str[-1]
        med = df.groupby("kinase", observed=True)["rank_in_site"].median().sort_values()
        rank_scoreable = set(med.index)

        universe = ksea_scoreable & rank_scoreable
        add("=" * 78)
        add("NETWORK: {} (spectral)".format(net))
        add("=" * 78)
        add("  kinases ranked              : {}".format(len(rank_scoreable)))
        add("  comparison universe         : {}".format(len(universe)))
        add("    (kinases BOTH methods could nominate)")
        ksea_u = ksea_sig & universe
        add("  KSEA-significant within it  : {}".format(len(ksea_u)))
        add("")

        for cut in CUTOFFS:
            top = set(med.head(cut).index) & universe
            inter = top & ksea_u
            union = top | ksea_u
            N = len(universe)
            K = len(ksea_u)
            n = len(top)
            x = len(inter)
            # hypergeometric: P(X <= x) is the disjointness-supporting tail
            p_lower = stats.hypergeom.cdf(x, N, K, n)
            p_upper = stats.hypergeom.sf(x - 1, N, K, n)
            expected = n * K / N if N else np.nan
            rows.append({
                "network": net, "cutoff": cut, "universe": N,
                "ksea_sig": K, "ranking_top": n, "overlap": x,
                "expected_overlap": expected,
                "jaccard": (x / len(union)) if union else np.nan,
                "p_depletion": p_lower, "p_enrichment": p_upper,
            })
            add("  --- ranking top {} by median rank ---".format(cut))
            add("    ranking set in universe : {}".format(n))
            add("    overlap with KSEA       : {}".format(x))
            add("    overlap expected by chance: {:.2f}".format(expected))
            add("    Jaccard                 : {:.3f}".format(
                (x / len(union)) if union else float("nan")))
            add("    hypergeometric p, depletion (supports disjointness) : {:.4f}".format(
                p_lower))
            add("    hypergeometric p, enrichment                        : {:.4f}".format(
                p_upper))
            if x:
                add("    shared: " + ", ".join(sorted(inter)))
            add("")

        add("  top 20 kinases by median rank ({}): ".format(net))
        add("    " + ", ".join(med.head(20).index.tolist()))
        add("")

    res = pd.DataFrame(rows)
    if not res.empty:
        res.to_csv(OUT / "ksea_ranking_overlap.csv", index=False)

    add("=" * 78)
    add("VERDICT")
    add("=" * 78)
    if not res.empty:
        worst = res["overlap"].max()
        add("  Maximum overlap across all networks and cutoffs: {} kinase(s).".format(
            int(worst)))
        sig_dep = res[res["p_depletion"] < 0.05]
        add("  Comparisons where the depletion is significant at p < 0.05: {} of {}.".format(
            len(sig_dep), len(res)))
        add("")
        if worst == 0:
            add("  The two sets are DISJOINT at every cutoff examined. The Chapter II")
            add("  sentence is supported as written.")
        elif res["overlap"].max() <= 2:
            add("  The sets overlap in at most {} kinase(s). 'Largely disjoint' is".format(
                int(worst)))
            add("  supportable, but the sentence should state the cutoff and the")
            add("  overlap count rather than asserting disjointness in general.")
        else:
            add("  The overlap is larger than 'largely disjoint' implies. The")
            add("  Chapter II sentence should be softened to what the table shows.")
    add("")
    add("Written: thesis/analysis/out/ksea_ranking_overlap.{txt,csv}")

    text = "\n".join(L)
    (OUT / "ksea_ranking_overlap.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
