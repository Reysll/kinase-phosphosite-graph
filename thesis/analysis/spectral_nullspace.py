"""
What the wasted spectral dimensions actually contain: sqrt(degree).

WHY
---
`spectral_spectrum.py` found that of the 33 eigenvalues
`SpectralEmbeddingStrategy` computes, the number that are numerically zero is
13 (generic), 18 (control) and 17 (liver FC). The first genuinely non-zero
eigenvalue does not appear until index 13, 18 and 17 respectively. The code
drops exactly ONE, so 12, 17 and 16 of the 32 shipped dimensions are
null-space vectors.

The normalized Laplacian L = I - D^-1/2 A D^-1/2 has null space spanned by the
component indicators scaled by sqrt(degree): for component C the vector
u_C(i) = sqrt(deg_i) if i in C else 0 satisfies L u_C = 0. ARPACK returns an
arbitrary orthonormal basis of that space, so each returned null vector is a
linear combination of the u_C. Restricted to the nodes of ONE component, every
such combination is a scalar multiple of sqrt(deg_i).

If that holds, then 16 of the 32 liver FC dimensions are, for the 98.5% of
nodes in the giant component, rescaled copies of sqrt(node degree). Three of
the project's central results would follow mechanically:

  - the degree-preserving null reproduces 96% of ranks (finding 9), because
    half the embedding is degree and the null preserves degree exactly;
  - hub bias persists under spectral (Phase 2), although the class docstring
    claims "High-degree nodes are penalized by construction ... Spectral
    embeddings do not share this property";
  - the site cloud collapses (rho = 0.9160), because sites in the giant
    component differ only in degree across those dimensions.

This script tests the claim directly and measures what a fix recovers. Two
quantities per network:

  1. R^2 of each shipped dimension regressed on sqrt(deg), over giant-component
     nodes only, and the share of sqrt(deg)'s norm that lies inside the span of
     the null dimensions.
  2. The counterfactual embedding: restrict the graph to its giant component,
     where the null space is one-dimensional, so dropping one eigenvector is
     correct and all 32 dimensions are informative. Reports the resulting site
     and kinase rho, and how many candidate kinases and FC-valid sites fall
     outside the giant component and would lose their embedding.

Read-only. Writes only into thesis/analysis/out/relational_only/.

Usage
-----
    python -m thesis.analysis.spectral_nullspace
"""

from __future__ import annotations

import json
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from thesis.analysis.embedding_geometry import _graph_for
from thesis.analysis.spectral_spectrum import ZERO_TOL, build_laplacian

OUT = Path("thesis/analysis/out/relational_only")


def _embedding_graph(network: str):
    from src_prediction.embeddings.embeddings import build_graph_from_edges
    from src_prediction.engine.leave_one_out import _build_embedding_relations
    graph, allowed = _graph_for(network)
    rels = _build_embedding_relations(graph.edges, allowed)
    eg = build_graph_from_edges(
        edges_df=graph.edges, allowed_relations=rels, directed=False)
    return eg.to_undirected() if eg.is_directed() else eg


def _rho(M):
    mu = M.mean(axis=0)
    dev = M - mu
    a = float(mu @ mu)
    b = float((dev * dev).sum(axis=1).mean())
    return float(np.sqrt(a / (a + b))) if a + b > 0 else float("nan")


def _spectral(L, n, dims, drop):
    """Compute `dims` eigenvectors after dropping the first `drop`."""
    k = min(dims + drop, n - 1)
    rng = np.random.default_rng(42)
    evals, evecs = spla.eigsh(L, k=k, which="SM", v0=rng.standard_normal(n))
    order = np.argsort(evals)
    return evals[order], evecs[:, order[drop:drop + dims]]


def analyse(network: str, site_ids, kinase_ids) -> dict:
    und = _embedding_graph(network)
    nodes, L, deg = build_laplacian(und)
    n = len(nodes)
    idx = {str(nd): i for i, nd in enumerate(nodes)}

    giant = max(nx.connected_components(und), key=len)
    giant_rows = np.array(sorted(idx[str(x)] for x in giant))

    evals, _ = _spectral(L, n, 32, 1)
    n_zero = int((evals[:33] < ZERO_TOL).sum())
    shipped = _spectral(L, n, 32, 1)[1]

    # --- 1. is the null space sqrt(degree) on the giant component?
    s = np.sqrt(deg[giant_rows])
    s = s - s.mean()
    s_n = s / np.linalg.norm(s)
    r2 = []
    for j in range(shipped.shape[1]):
        v = shipped[giant_rows, j]
        v = v - v.mean()
        nv = np.linalg.norm(v)
        r2.append(float((s_n @ (v / nv)) ** 2) if nv > 0 else float("nan"))

    n_null_dims = max(n_zero - 1, 0)
    if n_null_dims > 0:
        B = shipped[giant_rows][:, :n_null_dims]
        B = B - B.mean(axis=0)
        q, _ = np.linalg.qr(B)
        captured = float(np.linalg.norm(q.T @ s_n) ** 2)
    else:
        captured = float("nan")

    # --- 2. the fix: restrict to the giant component
    sub = und.subgraph(giant)
    sub_nodes, Lsub, _ = build_laplacian(sub)
    _, fixed = _spectral(Lsub, len(sub_nodes), 32, 1)
    fidx = {str(nd): i for i, nd in enumerate(sub_nodes)}

    srows_ship = [idx[s_] for s_ in site_ids if s_ in idx]
    krows_ship = [idx[k_] for k_ in kinase_ids if k_ in idx]
    srows_fix = [fidx[s_] for s_ in site_ids if s_ in fidx]
    krows_fix = [fidx[k_] for k_ in kinase_ids if k_ in fidx]

    row = {
        "network": network,
        "n_nodes": n,
        "giant_component": len(giant),
        "n_zero_eigenvalues": n_zero,
        "n_null_dims_shipped": n_null_dims,
        "n_informative_dims_shipped": 32 - n_null_dims,
        "r2_sqrtdeg_median_null_dims": float(np.median(r2[:n_null_dims])) if n_null_dims else float("nan"),
        "r2_sqrtdeg_max_null_dims": float(np.max(r2[:n_null_dims])) if n_null_dims else float("nan"),
        "r2_sqrtdeg_median_informative_dims": float(np.median(r2[n_null_dims:])) if n_null_dims < 32 else float("nan"),
        "sqrtdeg_norm_share_in_null_span": captured,
        "site_rho_shipped": _rho(shipped[srows_ship]),
        "kinase_rho_shipped": _rho(shipped[krows_ship]),
        "site_rho_giant_only": _rho(fixed[srows_fix]),
        "kinase_rho_giant_only": _rho(fixed[krows_fix]),
        "sites_in_graph": len(srows_ship),
        "sites_lost_by_giant_restriction": len(srows_ship) - len(srows_fix),
        "kinases_in_graph": len(krows_ship),
        "kinases_lost_by_giant_restriction": len(krows_ship) - len(krows_fix),
        "r2_per_dim": r2,
    }
    print("[null] {}: null dims {}/32  R2(sqrt deg) med {:.4f}  "
          "sqrtdeg in null span {:.4f}  site_rho {:.4f} -> {:.4f}".format(
              network, n_null_dims, row["r2_sqrtdeg_median_null_dims"],
              captured, row["site_rho_shipped"], row["site_rho_giant_only"]))
    return row


def main() -> None:
    from src_prediction.core.config import CANDIDATE_KINASES_OUT
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    OUT.mkdir(parents=True, exist_ok=True)
    kinase_ids = read_csv_auto(CANDIDATE_KINASES_OUT)["node_id"].astype(str).tolist()
    site_ids = sorted(load_site_fc_values(LIVER_EXCEL).keys())

    rows = [analyse(net, site_ids, kinase_ids)
            for net in ("generic", "control", "liver_fc")]

    pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, list)}
                  for r in rows]).to_csv(OUT / "spectral_nullspace.csv", index=False)
    (OUT / "spectral_nullspace_meta.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8")

    L = ["SPECTRAL NULL SPACE: the wasted dimensions are sqrt(degree)",
         "=" * 70, ""]
    for r in rows:
        L.append(r["network"].upper())
        L.append("  nodes {:,}  giant {:,}  zero eigenvalues {}".format(
            r["n_nodes"], r["giant_component"], r["n_zero_eigenvalues"]))
        L.append("  shipped dims: {} null, {} informative".format(
            r["n_null_dims_shipped"], r["n_informative_dims_shipped"]))
        L.append("  R^2 vs sqrt(deg) on giant component:")
        L.append("     null dims        median {:.4f}  max {:.4f}".format(
            r["r2_sqrtdeg_median_null_dims"], r["r2_sqrtdeg_max_null_dims"]))
        L.append("     informative dims median {:.4f}".format(
            r["r2_sqrtdeg_median_informative_dims"]))
        L.append("  share of sqrt(deg) inside the null-dim span: {:.4f}".format(
            r["sqrtdeg_norm_share_in_null_span"]))
        L.append("  FIX, restrict to giant component:")
        L.append("     site rho   {:.4f} -> {:.4f}".format(
            r["site_rho_shipped"], r["site_rho_giant_only"]))
        L.append("     kinase rho {:.4f} -> {:.4f}".format(
            r["kinase_rho_shipped"], r["kinase_rho_giant_only"]))
        L.append("     lost: {} of {} sites, {} of {} kinases".format(
            r["sites_lost_by_giant_restriction"], r["sites_in_graph"],
            r["kinases_lost_by_giant_restriction"], r["kinases_in_graph"]))
        L.append("")
    txt = "\n".join(L)
    (OUT / "spectral_nullspace.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
