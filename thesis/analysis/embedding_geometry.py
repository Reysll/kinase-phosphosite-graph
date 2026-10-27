r"""
Embedding geometry: can relational-only features (Option D) actually work?

WHY
---
Finding 26 is that the kinase ranking is almost perfectly site-invariant
(eta^2 = 1.0000 under uniform negatives). The leading hypothesis (STATE.md
open decision 2) is the K block in `[K, S, |K-S|, K*S]`: K is identical at
every site a kinase is scored against, so logistic regression can learn a
per-kinase constant and add it everywhere.

Option D drops the K and S blocks and keeps only `|K-S|` and `K*S`, so no
feature is literally constant across sites. But that does NOT guarantee the
per-kinase constant disappears, and the reason is geometric rather than
algorithmic. For a fixed kinase the Hadamard term contributes

    w . (K * S) = (w * K) . S

which is a linear functional of S with kinase-specific coefficients. Write
S = S_bar + d, where S_bar is the mean site embedding and d the site's
deviation from it. Then

    (w * K) . S = (w * K) . S_bar  +  (w * K) . d
                  \_____________/     \_________/
                   per-kinase CONSTANT   site-varying part

So a per-kinase base rate survives Option D, carried by the mean site
embedding, and the share of the score it controls is governed by how
concentrated the site embeddings are. The relevant quantity is

    rho = ||S_bar|| / sqrt(||S_bar||^2 + E||d||^2)

the fraction of site-embedding magnitude that is common to every site. If
rho is near 1 the sites are nearly identical as far as the model can see,
the constant part dominates, and Option D will NOT restore site dependence.
If rho is small the deviations dominate and Option D has room to work.

This is a read-only diagnostic. It builds the same embeddings the runs use
and measures the geometry. It costs minutes, not hours, and it makes Option D
a falsifiable prediction instead of a hope.

Writes only into thesis/analysis/out/relational_only/.

Usage
-----
    python -m thesis.analysis.embedding_geometry
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

OUT = Path("thesis/analysis/out/relational_only")


def _graph_for(network: str):
    from src_prediction.core.config import (
        GENERIC_EDGES, GENERIC_NODES, LIVER_EDGES, LIVER_NODES)
    from src_prediction.core.graph_loader import load_graph
    from src_prediction.core.relation_filters import (
        GENERIC_BASE_RELATIONS, SITE_CORR_CTRL_RELATIONS, SITE_CORR_FC_RELATIONS)

    if network == "generic":
        g = load_graph(GENERIC_NODES, GENERIC_EDGES, label="generic")
        return g, set(GENERIC_BASE_RELATIONS)
    g = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
    if network == "liver_fc":
        return g, set(GENERIC_BASE_RELATIONS) | set(SITE_CORR_FC_RELATIONS)
    if network == "control":
        return g, set(GENERIC_BASE_RELATIONS) | set(SITE_CORR_CTRL_RELATIONS)
    raise ValueError(network)


def embed(network: str):
    """Reproduce the embeddings the spectral runs use, PSP edges excluded."""
    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy
    from src_prediction.embeddings.embeddings import build_graph_from_edges
    from src_prediction.engine.leave_one_out import _build_embedding_relations

    graph, allowed = _graph_for(network)
    rels = _build_embedding_relations(graph.edges, allowed)
    strategy = SpectralEmbeddingStrategy(dimensions=32, directed=False, random_state=42)
    eg = build_graph_from_edges(
        edges_df=graph.edges, allowed_relations=rels, directed=strategy.directed)
    print("[geom] {}: embedding graph {:,} nodes, {:,} edges".format(
        network, eg.number_of_nodes(), eg.number_of_edges()))
    return strategy.fit(eg)


def geometry(emb: dict, kinase_ids, site_ids) -> dict:
    """Concentration of the site cloud and of the kinase cloud."""
    K = np.array([emb[k] for k in kinase_ids if k in emb])
    S = np.array([emb[s] for s in site_ids if s in emb])

    out = {"n_kinases": int(len(K)), "n_sites": int(len(S))}
    for name, M in (("site", S), ("kinase", K)):
        mu = M.mean(axis=0)
        dev = M - mu
        mu_sq = float(mu @ mu)
        dev_sq = float((dev * dev).sum(axis=1).mean())
        rho = float(np.sqrt(mu_sq / (mu_sq + dev_sq))) if mu_sq + dev_sq > 0 else float("nan")
        out["{}_mean_norm".format(name)] = float(np.linalg.norm(mu))
        out["{}_rms_deviation".format(name)] = float(np.sqrt(dev_sq))
        out["{}_rho_common".format(name)] = rho
    return out


def hadamard_split(emb: dict, kinase_ids, site_ids, n_sites=400, seed=0) -> dict:
    """
    How much of the Hadamard term is the per-kinase constant?

    For each kinase, compute K*S over a sample of sites and split the result
    into the part common to all sites (K*S_bar) and the site-varying part.
    Reported as the median over kinases of
        ||K*S_bar|| / sqrt(||K*S_bar||^2 + E||K*(S-S_bar)||^2)
    which is the share of the Hadamard feature vector's magnitude that no
    choice of logistic-regression weights can make site-dependent.
    """
    rng = np.random.default_rng(seed)
    sites = [s for s in site_ids if s in emb]
    if len(sites) > n_sites:
        sites = [sites[i] for i in rng.choice(len(sites), size=n_sites, replace=False)]
    S = np.array([emb[s] for s in sites])
    S_bar = S.mean(axis=0)
    D = S - S_bar

    shares = []
    for k in kinase_ids:
        if k not in emb:
            continue
        Kv = emb[k]
        const = Kv * S_bar
        const_sq = float(const @ const)
        var_sq = float(((Kv * D) ** 2).sum(axis=1).mean())
        if const_sq + var_sq > 0:
            shares.append(np.sqrt(const_sq / (const_sq + var_sq)))
    shares = np.array(shares)
    return {
        "hadamard_const_share_median": float(np.median(shares)),
        "hadamard_const_share_p10": float(np.percentile(shares, 10)),
        "hadamard_const_share_p90": float(np.percentile(shares, 90)),
        "n_sites_sampled": int(len(sites)),
        "n_kinases_scored": int(len(shares)),
    }


def main() -> None:
    from src_prediction.core.config import CANDIDATE_KINASES_OUT
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    OUT.mkdir(parents=True, exist_ok=True)
    cands = read_csv_auto(CANDIDATE_KINASES_OUT)["node_id"].astype(str).tolist()
    target_sites = sorted(load_site_fc_values(LIVER_EXCEL).keys())
    print("[geom] {:,} candidate kinases, {:,} FC-valid sites".format(
        len(cands), len(target_sites)))

    rows = []
    for network in ("generic", "control", "liver_fc"):
        emb = embed(network)
        row = {"network": network}
        row.update(geometry(emb, cands, target_sites))
        row.update(hadamard_split(emb, cands, target_sites))
        rows.append(row)
        print("[geom] {}: site_rho={:.4f} kinase_rho={:.4f} hadamard_const={:.4f}".format(
            network, row["site_rho_common"], row["kinase_rho_common"],
            row["hadamard_const_share_median"]))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "embedding_geometry.csv", index=False)

    lines = ["EMBEDDING GEOMETRY: can Option D restore site dependence?",
             "=" * 70, ""]
    lines.append("rho_common is the fraction of embedding magnitude shared by every")
    lines.append("node of that type. High site rho means the sites look alike to the")
    lines.append("model, so a per-kinase constant survives relational-only features.")
    lines.append("")
    for r in rows:
        lines.append("{}".format(r["network"]).upper())
        lines.append("  sites   : n={:,}  ||mean||={:.4f}  rms dev={:.4f}  rho={:.4f}".format(
            r["n_sites"], r["site_mean_norm"], r["site_rms_deviation"],
            r["site_rho_common"]))
        lines.append("  kinases : n={:,}  ||mean||={:.4f}  rms dev={:.4f}  rho={:.4f}".format(
            r["n_kinases"], r["kinase_mean_norm"], r["kinase_rms_deviation"],
            r["kinase_rho_common"]))
        lines.append("  Hadamard K*S constant share (median over kinases): {:.4f}".format(
            r["hadamard_const_share_median"]))
        lines.append("     p10 {:.4f}  p90 {:.4f}  over {:,} sites x {:,} kinases".format(
            r["hadamard_const_share_p10"], r["hadamard_const_share_p90"],
            r["n_sites_sampled"], r["n_kinases_scored"]))
        lines.append("")
    txt = "\n".join(lines)
    (OUT / "embedding_geometry.txt").write_text(txt, encoding="utf-8")
    (OUT / "embedding_geometry_meta.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
