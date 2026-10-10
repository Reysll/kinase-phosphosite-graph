"""
Tier 1 verification analyses (read-only).

T1.1  Within-category ranking: stratified (rich/poor) and kinase-level paired
      tests on adjusted_held_out_rank_in_category.
      SCOPE LIMIT: that column exists only in the four SPECTRAL runs (Phase 1
      node2vec scripts never passed kinase_categories), so the family here is
      the 3 spectral context tests. The embedding block cannot be run.
      Goal: separate "poor kinases lose their own signal" from "poor kinases
      are displaced in global rank by rich kinases".

T1.2  Context-network similarity: Spearman correlation of per-trial ranks
      between every network pair, plus top-10-per-site overlap from the
      Phase 3 predict-all matrices. Run for spectral and node2vec.

T1.3  Correlation-edge accounting: how many edges each context network adds,
      and how those edge endpoints distribute over rich/average/poor kinase
      substrates.

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.tier1_analyses

Seed: none required (deterministic).
"""

from __future__ import annotations

import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from thesis.analysis.stat_tests import KEY, holm, load_all, paired_test, union_unrankable
from thesis.analysis.verify_tests import categories, fmt

OUT = Path("thesis/analysis/out")
PRED = Path("outputs_prediction")
CAT_COL = "adjusted_held_out_rank_in_category"
RANK_COL = "adjusted_held_out_rank"
SPECTRAL_NETS = ["Generic", "Control", "Cancer", "Liver FC"]
DIRS = {"Generic": "generic", "Control": "control", "Cancer": "cancer", "Liver FC": "liver"}
INFER = {"Generic": "generic", "Control": "control", "Cancer": "cancer", "Liver FC": "liver_fc"}


# ------------------------------------------------------------------ T1.1
def load_category_frames() -> dict[str, pd.DataFrame]:
    out = {}
    for net, d in DIRS.items():
        p = PRED / f"{d}_multi_kinase_spectral" / "results.csv.gz"
        out[net] = pd.read_csv(p)[
            [KEY, "true_kinase_node_id", RANK_COL, CAT_COL, "held_out_kinase_category"]
        ]
    return out


def t11(keep: set[int], cats: dict[str, str], L: list[str]) -> None:
    add = L.append
    frames = load_category_frames()
    add("T1.1  WITHIN-CATEGORY RANKING (adjusted_held_out_rank_in_category)")
    add("=" * 78)
    add("SCOPE: spectral runs only -- Phase 1 node2vec results lack this column.")
    add("Family = 3 context tests (Generic vs Control/Cancer/Liver FC), Holm-corrected")
    add("within each stratum. Convention: positive diff / positive r = B better.")
    add("")

    for stratum in ("rich", "poor", "average"):
        ks = {k for k, c in cats.items() if c == stratum}
        for level in ("trial", "kinase"):
            rows = []
            for other in ("Control", "Cancer", "Liver FC"):
                da = frames["Generic"]
                db = frames[other]
                da = da[da[KEY].isin(keep) & da["true_kinase_node_id"].isin(ks)]
                db = db[db[KEY].isin(keep) & db["true_kinase_node_id"].isin(ks)]
                m = da.merge(db, on=KEY, suffixes=("_a", "_b"))
                m = m.dropna(subset=[f"{CAT_COL}_a", f"{CAT_COL}_b"])
                if m.empty:
                    continue
                if level == "kinase":
                    g = m.groupby("true_kinase_node_id_a")[
                        [f"{CAT_COL}_a", f"{CAT_COL}_b"]
                    ].median()
                    a, b = g[f"{CAT_COL}_a"].to_numpy(), g[f"{CAT_COL}_b"].to_numpy()
                else:
                    a, b = m[f"{CAT_COL}_a"].to_numpy(), m[f"{CAT_COL}_b"].to_numpy()
                res = paired_test(a, b, "Generic", other)
                d = a - b
                res["n_B_better"] = int((d > 0).sum())
                res["n_A_better"] = int((d < 0).sum())
                rows.append(res)
            if not rows:
                continue
            t = pd.DataFrame(rows)
            t["p_holm"] = holm(t["p_value"].tolist())
            t["significant_holm_0.05"] = t["p_holm"] < 0.05
            t = t[["comparison", "n", "n_ties", "median_A", "median_B",
                   "median_paired_diff", "rank_biserial", "n_B_better", "n_A_better",
                   "p_value", "p_holm", "significant_holm_0.05", "favours"]]
            add(f"--- {stratum.upper()} | {level} level | within-category rank ---")
            add(fmt(t))
            add("")
            t.to_csv(OUT / f"t11_withincat_{stratum}_{level}.csv", index=False)

    # side-by-side: global vs within-category median for the same trials
    add("Global vs within-category median rank (spectral, same trials)")
    add("-" * 78)
    rows = []
    for stratum in ("rich", "average", "poor"):
        ks = {k for k, c in cats.items() if c == stratum}
        r = {"category": stratum}
        for net in SPECTRAL_NETS:
            d = frames[net]
            d = d[d[KEY].isin(keep) & d["true_kinase_node_id"].isin(ks)]
            r[f"{net} global"] = float(d[RANK_COL].median())
            r[f"{net} within-cat"] = float(d[CAT_COL].median())
        rows.append(r)
    gv = pd.DataFrame(rows)
    add(gv.to_string(index=False))
    gv.to_csv(OUT / "t11_global_vs_withincat.csv", index=False)
    add("")


# ------------------------------------------------------------------ T1.2
def t12(frames, keep: set[int], L: list[str]) -> None:
    add = L.append
    add("T1.2  CONTEXT-NETWORK SIMILARITY")
    add("=" * 78)
    add("(a) Spearman rho of per-trial adjusted_held_out_rank between network pairs")
    add("")
    rows = []
    for emb in ("spectral", "node2vec"):
        for a, b in itertools.combinations(SPECTRAL_NETS, 2):
            da = frames[(a, emb)]
            db = frames[(b, emb)]
            m = (da[da[KEY].isin(keep)].merge(
                db[db[KEY].isin(keep)], on=KEY, suffixes=("_a", "_b")))
            rho, p = stats.spearmanr(m[f"{RANK_COL}_a"], m[f"{RANK_COL}_b"])
            ident = float((m[f"{RANK_COL}_a"] == m[f"{RANK_COL}_b"]).mean())
            rows.append({"embedding": emb, "pair": f"{a} vs {b}", "n": len(m),
                         "spearman_rho": rho, "p_value": p,
                         "frac_identical_rank": ident})
    sp = pd.DataFrame(rows)
    disp = sp.copy()
    disp["spearman_rho"] = disp["spearman_rho"].map("{:.4f}".format)
    disp["frac_identical_rank"] = disp["frac_identical_rank"].map("{:.3f}".format)
    disp["p_value"] = disp["p_value"].map(lambda v: f"{v:.2e}")
    add(disp.to_string(index=False))
    sp.to_csv(OUT / "t12_spearman.csv", index=False)
    add("")

    add("(b) Top-10-per-site overlap (Phase 3 predict-all, 3,228 FC-valid sites)")
    add("    mean |intersection| of each site's top-10 kinase set, out of 10")
    add("")
    rows = []
    for emb in ("spectral", "node2vec"):
        tops: dict[str, dict] = {}
        for net in SPECTRAL_NETS:
            p = PRED / "inference_all" / f"{INFER[net]}_{emb}" / "ranks.csv.gz"
            if not p.exists():
                continue
            df = pd.read_csv(p, usecols=["site_node_id", "kinase_node_id", "rank_in_site"])
            df = df[df["rank_in_site"] <= 10]
            tops[net] = {s: set(g) for s, g in
                         df.groupby("site_node_id")["kinase_node_id"]}
        for a, b in itertools.combinations([n for n in SPECTRAL_NETS if n in tops], 2):
            common = set(tops[a]) & set(tops[b])
            if not common:
                continue
            ov = [len(tops[a][s] & tops[b][s]) for s in common]
            rows.append({"embedding": emb, "pair": f"{a} vs {b}",
                         "n_sites": len(common), "mean_overlap_of_10": float(np.mean(ov)),
                         "frac_sites_identical_top10": float(np.mean([o == 10 for o in ov]))})
    if rows:
        tv = pd.DataFrame(rows)
        disp = tv.copy()
        disp["mean_overlap_of_10"] = disp["mean_overlap_of_10"].map("{:.2f}".format)
        disp["frac_sites_identical_top10"] = disp["frac_sites_identical_top10"].map("{:.3f}".format)
        add(disp.to_string(index=False))
        tv.to_csv(OUT / "t12_top10_overlap.csv", index=False)
    else:
        add("    (no inference_all matrices found)")
    add("")


# ------------------------------------------------------------------ T1.3
def t13(cats: dict[str, str], L: list[str]) -> None:
    add = L.append
    add("T1.3  CORRELATION-EDGE ACCOUNTING")
    add("=" * 78)
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.core.config import GENERIC_EDGES, LIVER_EDGES

    liver = read_csv_auto(LIVER_EDGES)
    counts = liver["relation"].value_counts()
    corr_sets = {
        "Control": ["site_corr_ctrl_pos", "site_corr_ctrl_neg"],
        "Cancer": ["site_corr_cancer_pos", "site_corr_cancer_neg"],
        "Liver FC": ["site_corr_fc_pos", "site_corr_fc_neg"],
    }
    base = int(counts[[r for r in counts.index if not r.startswith("site_corr")]].sum())
    add(f"Base (non-correlation) edges in liver file : {base:,}")
    rows = []
    for net, rels in corr_sets.items():
        pos, neg = int(counts.get(rels[0], 0)), int(counts.get(rels[1], 0))
        rows.append({"network": net, "pos_edges": pos, "neg_edges": neg,
                     "corr_edges_added": pos + neg,
                     "total_network_edges": base + pos + neg})
    et = pd.DataFrame(rows)
    add(et.to_string(index=False))
    et.to_csv(OUT / "t13_edge_counts.csv", index=False)
    add("")

    # site -> categories of its PSP kinases
    gen = read_csv_auto(GENERIC_EDGES)
    psp = gen[gen["relation"] == "phosphorylates"]
    site_cat: dict[str, set[str]] = {}
    for k, s in zip(psp["source"].astype(str), psp["target"].astype(str)):
        c = cats.get(k)
        if c:
            site_cat.setdefault(s, set()).add(c)

    add("Correlation-edge endpoints by substrate annotation")
    add("  'rich-touching'  = >=1 endpoint is a substrate of a rich kinase")
    add("  'poor-touching'  = >=1 endpoint is a substrate of a poor kinase")
    add("  'unannotated'    = NEITHER endpoint has any PSP kinase")
    add("")
    rows = []
    for net, rels in corr_sets.items():
        sub = liver[liver["relation"].isin(rels)]
        u = sub["source"].astype(str).to_numpy()
        v = sub["target"].astype(str).to_numpy()
        n = len(sub)
        def touch(cat: str) -> int:
            return int(sum(1 for a, b in zip(u, v)
                           if cat in site_cat.get(a, ()) or cat in site_cat.get(b, ())))
        unann = int(sum(1 for a, b in zip(u, v)
                        if not site_cat.get(a) and not site_cat.get(b)))
        both_ann = int(sum(1 for a, b in zip(u, v)
                           if site_cat.get(a) and site_cat.get(b)))
        rows.append({
            "network": net, "corr_edges": n,
            "rich_touching": touch("rich"), "rich_pct": 100 * touch("rich") / n,
            "average_touching": touch("average"),
            "poor_touching": touch("poor"), "poor_pct": 100 * touch("poor") / n,
            "both_endpoints_annotated": both_ann,
            "neither_annotated": unann, "neither_pct": 100 * unann / n,
        })
    ep = pd.DataFrame(rows)
    disp = ep.copy()
    for c in ("rich_pct", "poor_pct", "neither_pct"):
        disp[c] = disp[c].map("{:.1f}%".format)
    add(disp.to_string(index=False))
    ep.to_csv(OUT / "t13_endpoint_categories.csv", index=False)
    add("")
    n_sites_ann = len(site_cat)
    add(f"Sites with >=1 PSP kinase: {n_sites_ann:,}")
    for c in ("rich", "average", "poor"):
        add(f"  sites with a {c:<7} kinase: {sum(1 for v in site_cat.values() if c in v):,}")
    add("")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames, _ = load_all()
    keep = set(frames[("Generic", "node2vec")][KEY].astype(int)) - union_unrankable(frames)
    cats = categories()

    L: list[str] = ["TIER 1 ANALYSES", "=" * 78,
                    f"Trials retained: {len(keep):,}", ""]
    t11(keep, cats, L)
    t12(frames, keep, L)
    t13(cats, L)

    text = "\n".join(L)
    (OUT / "tier1_analyses.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
