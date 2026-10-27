"""
Does graph degree alone predict the kinase ranking?

WHY
---
`spectral_nullspace.py` established that a large share of the spectral
embedding is sqrt(node degree): on the giant component a null-space
eigenvector of L = I - D^-1/2 A D^-1/2 is a scalar multiple of sqrt(deg), the
context graphs have 116 and 122 components of size >= 2, and
`SpectralEmbeddingStrategy` drops only one of those null directions.

If the embedding is substantially degree, then three separate findings should
reduce to one mechanical cause:

  finding  9  the degree-preserving null reproduces 96% of ranks
  finding 13  eta^2 = 0.9997, the ranking is near-constant across sites
  finding 14  ranking kinases by PSP substrate count beats every model

This script tests that directly, against the Phase 3 output already on disk.
For each network it regresses a kinase's mean `rank_in_site` on its degree in
the embedding graph, and compares that against its PSP substrate count, which
is the popularity baseline's ranking key. Three questions:

  1. how much of the global ordering is explained by degree alone;
  2. how much by substrate count alone;
  3. whether degree explains the ordering better than substrate count does,
     which would make the baseline's strength a consequence of the embedding
     defect rather than an independent fact about the data.

Degree is computed on the embedding graph, which excludes `phosphorylates`, so
degree carries no direct information about the prediction target. Substrate
count is computed from the generic edge table, the same way
`kinase_categories` does it.

Read-only. Reads `outputs_prediction/inference_all/` and writes only into
thesis/analysis/out/relational_only/.

Usage
-----
    python -m thesis.analysis.degree_rank
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from thesis.analysis.spectral_nullspace import _embedding_graph

OUT = Path("thesis/analysis/out/relational_only")
PHASE3 = Path("outputs_prediction/inference_all")
NETWORKS = ("generic", "control", "liver_fc")
EMBEDDINGS = ("spectral", "node2vec", "concat")


def _ranks(network: str, embedding: str) -> pd.DataFrame | None:
    p = PHASE3 / "{}_{}".format(network, embedding) / "ranks.csv.gz"
    if not p.exists():
        print("[deg] missing {}".format(p))
        return None
    return pd.read_csv(p, usecols=["site_node_id", "kinase_node_id", "rank_in_site"])


def _fit(x: np.ndarray, y: np.ndarray) -> dict:
    ok = np.isfinite(x) & np.isfinite(y)
    x, y = x[ok], y[ok]
    if len(x) < 3:
        return {"n": int(len(x)), "pearson_r2": float("nan"),
                "spearman_rho": float("nan"), "spearman_p": float("nan")}
    r = stats.pearsonr(x, y)
    sp = stats.spearmanr(x, y)
    return {"n": int(len(x)), "pearson_r2": float(r.statistic ** 2),
            "spearman_rho": float(sp.statistic), "spearman_p": float(sp.pvalue)}


def analyse(network: str, embedding: str, deg: dict, sub_counts: dict) -> dict | None:
    ranks = _ranks(network, embedding)
    if ranks is None:
        return None

    mean_rank = ranks.groupby("kinase_node_id")["rank_in_site"].agg(["mean", "std", "count"])
    mean_rank = mean_rank.rename(columns={"mean": "mean_rank", "std": "sd_rank",
                                          "count": "n_sites"})
    mean_rank["degree"] = [deg.get(k, np.nan) for k in mean_rank.index]
    mean_rank["sqrt_degree"] = np.sqrt(mean_rank["degree"])
    mean_rank["psp_substrates"] = [sub_counts.get(k, 0) for k in mean_rank.index]
    mean_rank["log1p_substrates"] = np.log1p(mean_rank["psp_substrates"])
    mean_rank["network"] = network
    mean_rank["embedding"] = embedding
    mean_rank.to_csv(OUT / "degree_rank_{}_{}.csv".format(network, embedding))

    y = mean_rank["mean_rank"].values
    row = {"network": network, "embedding": embedding,
           "n_kinases": int(len(mean_rank)),
           "n_sites": int(ranks["site_node_id"].nunique())}
    for name, col in (("degree", "degree"), ("sqrt_degree", "sqrt_degree"),
                      ("substrates", "psp_substrates"),
                      ("log_substrates", "log1p_substrates")):
        f = _fit(mean_rank[col].values.astype(float), y)
        row["{}_r2".format(name)] = f["pearson_r2"]
        row["{}_rho".format(name)] = f["spearman_rho"]
        row["{}_p".format(name)] = f["spearman_p"]
    row["median_sd_rank"] = float(mean_rank["sd_rank"].median())
    print("[deg] {:<9} {:<9}: degree rho={:+.4f} R2={:.4f} | substrates "
          "rho={:+.4f} R2={:.4f}".format(
              network, embedding, row["degree_rho"], row["degree_r2"],
              row["substrates_rho"], row["substrates_r2"]))
    return row


def degree_composition(kinase_ids) -> pd.DataFrame:
    """
    What is a kinase's degree made of?

    If the ranking is degree (above) and 90.6% of generic edges are MSigDB
    pathway cliques (GRAPH_WALKTHROUGH section 3.3), then the ranking is
    largely a kinase's pathway-annotation breadth. This counts, per relation,
    how many edges touch a candidate kinase, on the embedding graph's relation
    set (`phosphorylates` excluded).
    """
    from src_prediction.core.config import LIVER_EDGES, LIVER_NODES
    from src_prediction.core.graph_loader import load_graph
    from src_prediction.core.relation_filters import PSP_KSA_RELATIONS

    edges = load_graph(LIVER_NODES, LIVER_EDGES, label="liver").edges
    edges = edges[~edges["relation"].astype(str).isin(set(PSP_KSA_RELATIONS))]
    kin = set(kinase_ids)
    touches = (edges["source"].astype(str).isin(kin)
               | edges["target"].astype(str).isin(kin))
    sub = edges[touches]
    comp = (sub.groupby(sub["relation"].astype(str)).size()
            .sort_values(ascending=False).rename("edges_touching_a_kinase"))
    out = comp.to_frame()
    out["share"] = out["edges_touching_a_kinase"] / out["edges_touching_a_kinase"].sum()
    return out


def main() -> None:
    from src_prediction.core.config import CANDIDATE_KINASES_OUT, GENERIC_EDGES
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.data_prep.kinase_categories import compute_psp_substrate_counts

    OUT.mkdir(parents=True, exist_ok=True)
    counts = compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES))
    counts = dict(counts)
    n_cands = len(read_csv_auto(CANDIDATE_KINASES_OUT))
    print("[deg] {} candidate kinases, {} with substrate counts".format(
        n_cands, len(counts)))

    rows = []
    for net in NETWORKS:
        deg = dict(_embedding_graph(net).degree())
        for emb in EMBEDDINGS:
            r = analyse(net, emb, deg, counts)
            if r:
                rows.append(r)
    if not rows:
        print("[deg] no Phase 3 output found; nothing to do")
        return

    pd.DataFrame(rows).to_csv(OUT / "degree_rank.csv", index=False)
    (OUT / "degree_rank_meta.json").write_text(json.dumps(rows, indent=2),
                                               encoding="utf-8")

    L = ["DEGREE AND RANK: is the global ordering just node degree?",
         "=" * 70, "",
         "Mean rank_in_site per kinase from the Phase 3 predict-all output,",
         "regressed on embedding-graph degree and on PSP substrate count.",
         "Degree is computed on the embedding graph, which excludes",
         "`phosphorylates`, so it carries no direct target information.", ""]
    hdr = "{:<9} {:<9} {:>8} {:>8} {:>9} {:>9} {:>8}".format(
        "network", "embedding", "deg_rho", "deg_R2", "subs_rho", "subs_R2", "sd_rank")
    L.append(hdr)
    L.append("-" * len(hdr))
    for r in rows:
        L.append("{:<9} {:<9} {:>+8.4f} {:>8.4f} {:>+9.4f} {:>9.4f} {:>8.2f}".format(
            r["network"], r["embedding"], r["degree_rho"], r["degree_r2"],
            r["substrates_rho"], r["substrates_r2"], r["median_sd_rank"]))
    L.append("")
    L.append("sd_rank is the median per-kinase SD of rank_in_site across sites;")
    L.append("small values mean the kinase holds the same rank at every site.")
    L.append("")

    comp = degree_composition(
        read_csv_auto(CANDIDATE_KINASES_OUT)["node_id"].astype(str).tolist())
    comp.to_csv(OUT / "degree_composition.csv")
    L.append("WHAT DEGREE IS MADE OF (liver edge table, PSP edges excluded,")
    L.append("edges touching a candidate kinase):")
    L.append("")
    L.append("  {:<28} {:>12} {:>8}".format("relation", "edges", "share"))
    for rel, r in comp.iterrows():
        L.append("  {:<28} {:>12,} {:>7.1%}".format(
            rel, int(r["edges_touching_a_kinase"]), r["share"]))
    L.append("")
    txt = "\n".join(L)
    (OUT / "degree_rank.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
