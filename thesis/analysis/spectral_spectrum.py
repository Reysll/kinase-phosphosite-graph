"""
Why the site embeddings collapse: eigenvalue multiplicity at zero.

WHY
---
`embedding_geometry.py` found that under Liver FC, 92% of the site-embedding
magnitude is common to all 3,228 sites (rho = 0.9160), and 95% of the Hadamard
feature vector is a per-kinase constant no weight vector can make
site-dependent. That is a sufficient explanation for eta^2 = 1.0000 on its own,
and it also means Option D cannot restore site dependence.

This script asks why the site cloud collapses. The suspect is structural, not
statistical. `SpectralEmbeddingStrategy.fit` does:

    k = min(self.dimensions + 1, n - 1)
    eigenvalues, eigenvectors = spla.eigsh(L, k=k, which="SM", ...)
    order = np.argsort(eigenvalues)[1 : self.dimensions + 1]

It requests 33 eigenvectors and drops exactly ONE, on the assumption that the
graph is connected and so has a single eigenvalue at zero. The normalized
Laplacian of a graph with c connected components has eigenvalue 0 with
multiplicity c. Every one of those c eigenvectors is an indicator that is
constant within each component. Drop only the first and the remaining c-1 enter
the embedding as near-constant columns, which inflates the common component of
every node in the giant component and wastes dimensions that should have
carried site structure.

Isolated nodes (degree 0) do not contribute here: `fit` sets d_inv_sqrt = 1 for
them, so their Laplacian row is the identity row and their eigenvalue is 1, not
0. Components of size >= 2 are what matter.

Reported, read-only:
  1. component structure of each embedding graph
  2. the 33 computed eigenvalues, and how many sit at zero
  3. per-dimension standard deviation across sites, which shows directly how
     many of the 32 columns are near-constant
  4. the counterfactual: site rho after discarding every zero eigenvector
     instead of just one, which is the cheapest available fix

Writes only into thesis/analysis/out/relational_only/.

Usage
-----
    python -m thesis.analysis.spectral_spectrum
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from thesis.analysis.embedding_geometry import _graph_for

OUT = Path("thesis/analysis/out/relational_only")
ZERO_TOL = 1e-8


def build_laplacian(eg):
    """Byte-for-byte the matrix SpectralEmbeddingStrategy.fit builds."""
    nodes = list(eg.nodes())
    n = len(nodes)
    node_to_idx = {nd: i for i, nd in enumerate(nodes)}
    src, dst = [], []
    for u, v in eg.edges():
        i, j = node_to_idx[u], node_to_idx[v]
        src.extend([i, j])
        dst.extend([j, i])
    A = sp.csr_matrix(
        (np.ones(len(src), dtype=np.float32), (src, dst)), shape=(n, n))
    A = (A + A.T).sign().astype(np.float32)
    deg = np.asarray(A.sum(axis=1)).flatten()
    d_inv_sqrt = sp.diags(np.where(deg > 0, deg ** -0.5, 1.0))
    L = sp.eye(n, format="csr") - d_inv_sqrt @ A @ d_inv_sqrt
    return nodes, L, deg


def analyse(network: str, site_ids, kinase_ids) -> dict:
    import networkx as nx
    from src_prediction.embeddings.embeddings import build_graph_from_edges
    from src_prediction.engine.leave_one_out import _build_embedding_relations

    graph, allowed = _graph_for(network)
    rels = _build_embedding_relations(graph.edges, allowed)
    eg = build_graph_from_edges(
        edges_df=graph.edges, allowed_relations=rels, directed=False)

    und = eg.to_undirected() if eg.is_directed() else eg
    comps = sorted((len(c) for c in nx.connected_components(und)), reverse=True)
    n_iso = sum(1 for c in comps if c == 1)
    comps_ge2 = [c for c in comps if c >= 2]

    nodes, L, deg = build_laplacian(eg)
    n = len(nodes)
    k = min(32 + 1, n - 1)
    rng = np.random.default_rng(42)
    evals, evecs = spla.eigsh(L, k=k, which="SM", v0=rng.standard_normal(n))
    order_all = np.argsort(evals)
    evals_sorted = evals[order_all]
    n_zero = int((evals_sorted < ZERO_TOL).sum())

    idx = {str(nd): i for i, nd in enumerate(nodes)}
    site_rows = [idx[s] for s in site_ids if s in idx]
    kin_rows = [idx[k_] for k_ in kinase_ids if k_ in idx]

    def rho(M):
        mu = M.mean(axis=0)
        dev = M - mu
        a = float(mu @ mu)
        b = float((dev * dev).sum(axis=1).mean())
        return float(np.sqrt(a / (a + b))) if a + b > 0 else float("nan")

    # as shipped: drop exactly one
    shipped = evecs[:, order_all[1:33]]
    # counterfactual: drop every zero eigenvector
    avail = order_all[n_zero:]
    fixed = evecs[:, avail] if len(avail) else evecs[:, order_all[1:33]]

    sd_sites = shipped[site_rows].std(axis=0)
    row = {
        "network": network,
        "n_nodes": n,
        "n_components_total": len(comps),
        "n_components_ge2": len(comps_ge2),
        "n_isolated_nodes": n_iso,
        "giant_component": comps[0] if comps else 0,
        "second_component": comps_ge2[1] if len(comps_ge2) > 1 else 0,
        "n_zero_eigenvalues": n_zero,
        "eig_min": float(evals_sorted[0]),
        "eig_1": float(evals_sorted[1]),
        "eig_2": float(evals_sorted[2]) if len(evals_sorted) > 2 else float("nan"),
        "eig_max_computed": float(evals_sorted[-1]),
        "n_dims_near_constant_sites": int((sd_sites < 1e-6).sum()),
        "site_rho_shipped": rho(shipped[site_rows]),
        "kinase_rho_shipped": rho(shipped[kin_rows]),
        "n_dims_after_dropping_all_zeros": int(fixed.shape[1]),
        "site_rho_dropping_all_zeros": rho(fixed[site_rows]),
        "kinase_rho_dropping_all_zeros": rho(fixed[kin_rows]),
        "eigenvalues": [float(x) for x in evals_sorted],
        "site_dim_sd": [float(x) for x in sd_sites],
        "components_ge2_top10": comps_ge2[:10],
    }
    print("[spec] {}: n={:,} comps(>=2)={} iso={:,} giant={:,} zeros={} "
          "site_rho {:.4f} -> {:.4f}".format(
              network, n, row["n_components_ge2"], n_iso, row["giant_component"],
              n_zero, row["site_rho_shipped"], row["site_rho_dropping_all_zeros"]))
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

    flat = pd.DataFrame([{k: v for k, v in r.items()
                          if not isinstance(v, list)} for r in rows])
    flat.to_csv(OUT / "spectral_spectrum.csv", index=False)
    (OUT / "spectral_spectrum_meta.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8")

    L = ["SPECTRAL SPECTRUM: why the site embeddings collapse", "=" * 70, ""]
    L.append("SpectralEmbeddingStrategy requests 33 eigenvectors and drops ONE,")
    L.append("assuming a single zero eigenvalue. A graph with c components of")
    L.append("size >= 2 has c of them, each constant within its component.")
    L.append("")
    for r in rows:
        L.append(r["network"].upper())
        L.append("  nodes {:,}   components(>=2) {}   isolated {:,}".format(
            r["n_nodes"], r["n_components_ge2"], r["n_isolated_nodes"]))
        L.append("  giant component {:,}   next largest {:,}".format(
            r["giant_component"], r["second_component"]))
        L.append("  zero eigenvalues among the 33 computed: {}".format(
            r["n_zero_eigenvalues"]))
        L.append("  smallest eigenvalues: {:.3e}  {:.3e}  {:.3e}".format(
            r["eig_min"], r["eig_1"], r["eig_2"]))
        L.append("  of the 32 shipped dims, near-constant across sites: {}".format(
            r["n_dims_near_constant_sites"]))
        L.append("  site rho    as shipped {:.4f}   dropping all zeros {:.4f} ({} dims)".format(
            r["site_rho_shipped"], r["site_rho_dropping_all_zeros"],
            r["n_dims_after_dropping_all_zeros"]))
        L.append("  kinase rho  as shipped {:.4f}   dropping all zeros {:.4f}".format(
            r["kinase_rho_shipped"], r["kinase_rho_dropping_all_zeros"]))
        L.append("")
    txt = "\n".join(L)
    (OUT / "spectral_spectrum.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
