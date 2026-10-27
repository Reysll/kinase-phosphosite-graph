"""
Which correlation-edge filter makes the site embedding informative?

WHY
---
Two separate defects were established upstream of this script.

  A. Null-space contamination. `SpectralEmbeddingStrategy` requests
     dimensions+1 eigenvectors and drops one, which is correct only for a
     connected graph. The embedding graphs have 34, 116 and 122 components of
     size >= 2, so 12, 17 and 16 of the 32 shipped dimensions are null-space
     vectors, each a scalar multiple of sqrt(deg) on the giant component
     (R^2 = 1.0000). See `spectral_fix.py`.

  B. Spectral-gap collapse. Even after removing the whole null space, the site
     cloud stays collapsed under the context networks: Liver FC site rho is
     0.8882 against 0.1239 for generic. The first non-zero eigenvalue is
     3.68e-03 for Liver FC against 3.53e-02 for generic, an order of magnitude
     smaller, so the context graph has no bottleneck for the embedding to
     resolve. With 92% of site-embedding magnitude common to all 3,228 sites,
     no feature construction downstream can be site-specific. This is the
     defect that actually blocks context-aware prediction.

Defect B points at how correlation edges are built: a pooled 80th-percentile
threshold on pairwise-complete Pearson correlations with `min_periods=6`, no
significance control, which adds ~976,855 edges to Liver FC, of which 77-79%
touch no annotated site (finding 8) and only 8% use all 18 samples. The
hypothesis is that these edges flood the graph and destroy its bottleneck
structure.

This script tests that hypothesis on the embedding alone, with no
leave-one-out run. For each candidate filter it rebuilds the embedding graph,
removes the null space properly, and measures whether sites become
distinguishable. The quantities that matter:

  - `first_nonzero_eigenvalue`, the spectral gap; larger means real structure
  - `site_rho_corrected`, the share of site-embedding magnitude common to
    every site; lower means sites are distinguishable
  - `degree_dims_corrected`, which should be 0 once the null space is gone
  - component count and site coverage, the cost of a stricter filter

Variants: the current 80th-percentile set, sample-count floors, a stricter
magnitude percentile, BH-FDR on the correlation p-value, and base relations
only. CLAUDE.md warns that base-relations-only on the liver graph "breaks
SpectralEmbeddingStrategy's eigensolver (near-zero eigenvalue multiplicity)";
that is defect A, which this script's eigensolver handles, so the variant is
measurable here.

Read-only: reads `outputs/` and writes only into
thesis/analysis/out/relational_only/. No project code or data is modified.

Usage
-----
    python -m thesis.analysis.corr_edge_variants [--k 96]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse.linalg as spla
from scipy import stats

from thesis.analysis.spectral_fix import TOL, DIMS, r2_vs
from thesis.analysis.spectral_nullspace import _rho
from thesis.analysis.spectral_spectrum import build_laplacian

OUT = Path("thesis/analysis/out/relational_only")
FC_RELATIONS = ("site_corr_fc_pos", "site_corr_fc_neg")


def corr_pvalue(r: np.ndarray, n: np.ndarray) -> np.ndarray:
    """Two-sided p for a Pearson r from n paired observations."""
    r = np.clip(np.abs(r), 0, 1 - 1e-12)
    df = n - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = r * np.sqrt(df / (1 - r ** 2))
    p = 2 * stats.t.sf(np.abs(t), df)
    return np.where(df > 0, p, np.nan)


def bh_fdr(p: np.ndarray) -> np.ndarray:
    ok = np.isfinite(p)
    out = np.full(len(p), np.nan)
    q = p[ok]
    order = np.argsort(q)
    ranked = q[order]
    m = len(ranked)
    adj = ranked * m / np.arange(1, m + 1)
    adj = np.minimum.accumulate(adj[::-1])[::-1]
    res = np.empty(m)
    res[order] = np.clip(adj, 0, 1)
    out[ok] = res
    return out


def variants(corr: pd.DataFrame) -> dict:
    """Boolean masks over the FC correlation edges."""
    w = corr["weight"].values.astype(float)
    n = corr["n_common"].values.astype(float)
    aw = np.abs(w)
    p = corr_pvalue(w, n)
    q = bh_fdr(p)
    thr_top25 = float(np.nanpercentile(aw, 75.0))  # top 25% of this pool = ~95th pct originally
    return {
        "fc_all_current": np.ones(len(corr), dtype=bool),
        "fc_n_common_ge12": n >= 12,
        "fc_n_common_eq18": n >= 18,
        "fc_abs_r_top25pct": aw >= thr_top25,
        "fc_fdr05": np.isfinite(q) & (q < 0.05),
        "fc_fdr05_and_n12": np.isfinite(q) & (q < 0.05) & (n >= 12),
        "base_only": np.zeros(len(corr), dtype=bool),
    }


def embed_and_measure(edges: pd.DataFrame, site_ids, kinase_ids, k: int) -> dict:
    """
    Measure the embedding of one edge-filter variant.

    A stricter filter fragments the graph, which inflates the Laplacian's
    zero-eigenvalue multiplicity. If that multiplicity exceeds the eigenvector
    budget there is nothing informative left to keep, so the embedding is
    recomputed on the giant component alone, where the null space is
    one-dimensional and `k = dims + 1` is sufficient. `scope` records which
    path was taken, since the two are not interchangeable: the giant-component
    embedding simply has no coordinates for nodes outside it.
    """
    from src_prediction.embeddings.embeddings import build_graph_from_edges
    rels = set(edges["relation"].astype(str).unique())
    eg = build_graph_from_edges(edges_df=edges, allowed_relations=rels, directed=False)
    und = eg.to_undirected() if eg.is_directed() else eg

    comps = [len(c) for c in nx.connected_components(und)]
    giant_nodes = max(nx.connected_components(und), key=len)
    n_zero_est = sum(1 for c in comps if c >= 2)
    scope = "full"
    if n_zero_est + DIMS > k:
        scope = "giant"
        und = und.subgraph(giant_nodes).copy()
        comps = [len(giant_nodes)]
        k = DIMS + 1
    nodes, L, deg = build_laplacian(und)
    n = len(nodes)
    idx = {str(nd): i for i, nd in enumerate(nodes)}
    grows = np.array(sorted(idx[str(x)] for x in giant_nodes))

    rng = np.random.default_rng(42)
    evals, evecs = spla.eigsh(L, k=min(k, n - 1), which="SM",
                              v0=rng.standard_normal(n))
    order = np.argsort(evals)
    evals, evecs = evals[order], evecs[:, order]
    n_zero = int((evals < TOL).sum())
    keep = np.where(evals >= TOL)[0][:DIMS]
    corrected = evecs[:, keep]
    shipped = evecs[:, 1:DIMS + 1]

    s = np.sqrt(deg[grows])
    s = s - s.mean()
    s_n = s / np.linalg.norm(s)
    srows = [idx[x] for x in site_ids if x in idx]
    krows = [idx[x] for x in kinase_ids if x in idx]

    return {
        "scope": scope,
        "n_nodes": n,
        "n_edges": int(und.number_of_edges()),
        "n_components_ge2": int(n_zero_est),
        "giant_component": int(len(giant_nodes)),
        "n_zero_eigenvalues": n_zero,
        "enough_eigs": bool(n_zero + DIMS <= len(evals)),
        "first_nonzero_eigenvalue": float(evals[n_zero]) if n_zero < len(evals) else float("nan"),
        "n_corrected_dims": int(corrected.shape[1]),
        "degree_dims_shipped": int(sum(
            1 for j in range(shipped.shape[1]) if r2_vs(shipped[grows, j], s_n) > 0.99)),
        "degree_dims_corrected": int(sum(
            1 for j in range(corrected.shape[1]) if r2_vs(corrected[grows, j], s_n) > 0.99)),
        "sites_covered": int(len(srows)),
        "kinases_covered": int(len(krows)),
        "site_rho_shipped": _rho(shipped[srows]) if srows else float("nan"),
        "site_rho_corrected": _rho(corrected[srows]) if srows else float("nan"),
        "kinase_rho_corrected": _rho(corrected[krows]) if krows else float("nan"),
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=96)
    ap.add_argument("--only", default=None,
                    help="comma-separated variant names; default is all")
    a = ap.parse_args()
    wanted = set(a.only.split(",")) if a.only else None

    from src_prediction.core.config import (
        CANDIDATE_KINASES_OUT, LIVER_EDGES, LIVER_NODES)
    from src_prediction.core.graph_loader import load_graph
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.core.relation_filters import (
        GENERIC_BASE_RELATIONS, PSP_KSA_RELATIONS)
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    OUT.mkdir(parents=True, exist_ok=True)
    kinase_ids = read_csv_auto(CANDIDATE_KINASES_OUT)["node_id"].astype(str).tolist()
    site_ids = sorted(load_site_fc_values(LIVER_EXCEL).keys())

    g = load_graph(LIVER_NODES, LIVER_EDGES, label="liver")
    edges = g.edges
    base_rels = set(GENERIC_BASE_RELATIONS) - set(PSP_KSA_RELATIONS)
    base = edges[edges["relation"].astype(str).isin(base_rels)]
    corr = edges[edges["relation"].astype(str).isin(FC_RELATIONS)].reset_index(drop=True)
    print("[corr] base {:,} edges, FC correlation {:,} edges".format(
        len(base), len(corr)))
    print("[corr] n_common: min {:.0f} median {:.0f} max {:.0f}; "
          "share at 18: {:.1%}".format(
              corr["n_common"].min(), corr["n_common"].median(),
              corr["n_common"].max(), (corr["n_common"] >= 18).mean()))

    rows = []
    for name, mask in variants(corr).items():
        if wanted and name not in wanted:
            continue
        kept = corr[mask]
        sub = pd.concat([base, kept], ignore_index=True)
        row = {"variant": name, "corr_edges_kept": int(len(kept)),
               "corr_edges_share": float(len(kept) / len(corr))}
        row.update(embed_and_measure(sub, site_ids, kinase_ids, a.k))
        rows.append(row)
        print("[corr] {:<20} corr {:>9,} ({:5.1%})  gap {:.3e}  "
              "site_rho {:.4f}  sites {:,}".format(
                  name, row["corr_edges_kept"], row["corr_edges_share"],
                  row["first_nonzero_eigenvalue"], row["site_rho_corrected"],
                  row["sites_covered"]))

    stem = "corr_edge_variants" if not wanted else "corr_edge_variants_extra"
    df = pd.DataFrame(rows)
    df.to_csv(OUT / "{}.csv".format(stem), index=False)
    (OUT / "{}_meta.json".format(stem)).write_text(
        json.dumps(rows, indent=2), encoding="utf-8")

    L = ["CORRELATION-EDGE VARIANTS: what restores a site-specific embedding",
         "=" * 78, "",
         "All variants use the same base relations plus a filtered subset of the",
         "Liver FC correlation edges. The embedding drops the whole null space,",
         "so `degree dims` should be 0 and differences are structural.",
         "",
         "gap   = first non-zero eigenvalue; larger means real bottleneck structure",
         "s_rho = share of site-embedding magnitude common to ALL sites;",
         "        lower means sites are distinguishable (generic reference: 0.1239)",
         ""]
    hdr = "{:<20} {:>9} {:>6} {:>5} {:>6} {:>10} {:>7} {:>7} {:>6}".format(
        "variant", "corr", "share", "comp", "scope", "gap", "s_rho", "k_rho", "sites")
    L.append(hdr)
    L.append("-" * len(hdr))
    for r in rows:
        L.append("{:<20} {:>9,} {:>5.1%} {:>5} {:>6} {:>10.3e} {:>7.4f} {:>7.4f} "
                 "{:>6,}".format(
                     r["variant"], r["corr_edges_kept"], r["corr_edges_share"],
                     r["n_components_ge2"], r["scope"],
                     r["first_nonzero_eigenvalue"], r["site_rho_corrected"],
                     r["kinase_rho_corrected"], r["sites_covered"]))
    L.append("")
    L.append("scope 'giant' means fragmentation exceeded the eigenvector budget")
    L.append("and the embedding was computed on the giant component alone, so")
    L.append("nodes outside it have no coordinates. Not interchangeable with 'full'.")
    L.append("")
    for r in rows:
        if not r["enough_eigs"]:
            L.append("WARNING {}: k too small, only {} corrected dims".format(
                r["variant"], r["n_corrected_dims"]))
        if r["degree_dims_corrected"]:
            L.append("WARNING {}: {} corrected dims still sqrt(deg)".format(
                r["variant"], r["degree_dims_corrected"]))
    txt = "\n".join(L)
    (OUT / "corr_edge_variants.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
