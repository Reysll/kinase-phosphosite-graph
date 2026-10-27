"""
How many spectral dimensions are sqrt(degree), and what dropping them recovers.

WHY
---
`spectral_nullspace.py` confirmed the mechanism: on the giant component, a
null-space eigenvector of L = I - D^-1/2 A D^-1/2 is a scalar multiple of
sqrt(deg), and the shipped generic embedding contains one such dimension with
R^2 = 1.0000 against sqrt(deg) while its other dimensions sit at R^2 = 0.0011.
For control, the MEDIAN of the 32 shipped dimensions has R^2 = 1.0000, so at
least half of that embedding is degree.

Two things were wrong with the earlier counting and are fixed here.

  1. The tolerance. `spectral_spectrum.py` used 1e-8, but the numerically-zero
     eigenvalues sit at 1.2e-08 to 1.8e-08, just above it, while the first
     genuinely non-zero eigenvalue is 3.5e-02 (generic) or 3.7e-03 (context).
     The gap spans five orders of magnitude, so 1e-6 separates them cleanly.

  2. Counting by eigenvalue is indirect. What matters is how many shipped
     dimensions are degree, so this script counts dimensions with R^2 > 0.99
     against sqrt(deg) over giant-component nodes directly.

It then builds the corrected embedding. `SpectralEmbeddingStrategy` requests
dimensions+1 eigenvectors and drops one, which is only right for a connected
graph. These graphs have 34, 116 and 122 components of size >= 2, so the null
space has that dimension and the shipped embedding is mostly degree. The fix
requests enough eigenvectors to find 32 with eigenvalue above tolerance, and
keeps only those.

Reported per network, as-shipped against corrected:
  - dimensions that are sqrt(deg) (R^2 > 0.99)
  - site and kinase rho, the share of embedding magnitude common to all nodes
    of that type
  - R^2 distribution against sqrt(deg)

This measures the embedding only. Whether the corrected embedding changes rank
predictions is a separate run and is not claimed here.

Read-only. Writes only into thesis/analysis/out/relational_only/.

Usage
-----
    python -m thesis.analysis.spectral_fix [--k 72]
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
import scipy.sparse.linalg as spla

from thesis.analysis.embedding_geometry import _graph_for
from thesis.analysis.spectral_nullspace import _embedding_graph, _rho
from thesis.analysis.spectral_spectrum import build_laplacian

OUT = Path("thesis/analysis/out/relational_only")
TOL = 1e-6
DIMS = 32


def r2_vs(v: np.ndarray, target_n: np.ndarray) -> float:
    v = v - v.mean()
    nv = np.linalg.norm(v)
    return float((target_n @ (v / nv)) ** 2) if nv > 0 else float("nan")


def analyse(network: str, site_ids, kinase_ids, k: int) -> dict:
    und = _embedding_graph(network)
    nodes, L, deg = build_laplacian(und)
    n = len(nodes)
    idx = {str(nd): i for i, nd in enumerate(nodes)}
    giant = max(nx.connected_components(und), key=len)
    grows = np.array(sorted(idx[str(x)] for x in giant))

    s = np.sqrt(deg[grows])
    s = s - s.mean()
    s_n = s / np.linalg.norm(s)

    t0 = time.time()
    rng = np.random.default_rng(42)
    evals, evecs = spla.eigsh(L, k=min(k, n - 1), which="SM",
                              v0=rng.standard_normal(n))
    order = np.argsort(evals)
    evals, evecs = evals[order], evecs[:, order]
    print("[fix] {}: {} eigenvalues in {:.1f} min".format(
        network, len(evals), (time.time() - t0) / 60.0))

    n_zero = int((evals < TOL).sum())
    shipped = evecs[:, 1:DIMS + 1]
    keep = np.where(evals >= TOL)[0][:DIMS]
    corrected = evecs[:, keep]

    r2_ship = [r2_vs(shipped[grows, j], s_n) for j in range(shipped.shape[1])]
    r2_corr = [r2_vs(corrected[grows, j], s_n) for j in range(corrected.shape[1])]

    srows = [idx[x] for x in site_ids if x in idx]
    krows = [idx[x] for x in kinase_ids if x in idx]

    row = {
        "network": network,
        "n_nodes": n,
        "n_components_ge2": sum(1 for c in nx.connected_components(und) if len(c) >= 2),
        "giant_component": len(giant),
        "n_eigenvalues_computed": int(len(evals)),
        "n_zero_eigenvalues": n_zero,
        "first_nonzero_eigenvalue": float(evals[n_zero]) if n_zero < len(evals) else float("nan"),
        "enough_computed": bool(n_zero + DIMS <= len(evals)),
        "n_corrected_dims": int(corrected.shape[1]),
        "shipped_dims_are_degree": int(sum(1 for x in r2_ship if x > 0.99)),
        "corrected_dims_are_degree": int(sum(1 for x in r2_corr if x > 0.99)),
        "shipped_r2_median": float(np.median(r2_ship)),
        "corrected_r2_median": float(np.median(r2_corr)),
        "site_rho_shipped": _rho(shipped[srows]),
        "site_rho_corrected": _rho(corrected[srows]),
        "kinase_rho_shipped": _rho(shipped[krows]),
        "kinase_rho_corrected": _rho(corrected[krows]),
        "r2_shipped": r2_ship,
        "r2_corrected": r2_corr,
        "eigenvalues": [float(x) for x in evals],
    }
    print("[fix] {}: zeros {}  degree dims {}/32 -> {}/{}  "
          "site_rho {:.4f} -> {:.4f}  kinase_rho {:.4f} -> {:.4f}".format(
              network, n_zero, row["shipped_dims_are_degree"],
              row["corrected_dims_are_degree"], row["n_corrected_dims"],
              row["site_rho_shipped"], row["site_rho_corrected"],
              row["kinase_rho_shipped"], row["kinase_rho_corrected"]))
    return row


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--k", type=int, default=72,
                    help="eigenvectors to compute; must exceed n_zero + 32")
    a = ap.parse_args()

    from src_prediction.core.config import CANDIDATE_KINASES_OUT
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.runners.run_inference_all import LIVER_EXCEL, load_site_fc_values

    OUT.mkdir(parents=True, exist_ok=True)
    kinase_ids = read_csv_auto(CANDIDATE_KINASES_OUT)["node_id"].astype(str).tolist()
    site_ids = sorted(load_site_fc_values(LIVER_EXCEL).keys())

    rows = [analyse(net, site_ids, kinase_ids, a.k)
            for net in ("generic", "control", "liver_fc")]

    pd.DataFrame([{k: v for k, v in r.items() if not isinstance(v, list)}
                  for r in rows]).to_csv(OUT / "spectral_fix.csv", index=False)
    (OUT / "spectral_fix_meta.json").write_text(
        json.dumps(rows, indent=2), encoding="utf-8")

    L = ["SPECTRAL FIX: dropping the null space instead of one eigenvector",
         "=" * 70, "",
         "Shipped: eigsh(k=33), drop 1, keep 32. Correct for a connected graph.",
         "Corrected: keep the first 32 eigenvectors with eigenvalue > {:g}.".format(TOL),
         "'degree dims' counts dimensions with R^2 > 0.99 against sqrt(deg)",
         "over giant-component nodes.", ""]
    for r in rows:
        L.append(r["network"].upper())
        L.append("  nodes {:,}  components(>=2) {}  giant {:,}".format(
            r["n_nodes"], r["n_components_ge2"], r["giant_component"]))
        L.append("  zero eigenvalues {} of {} computed; first non-zero {:.3e}".format(
            r["n_zero_eigenvalues"], r["n_eigenvalues_computed"],
            r["first_nonzero_eigenvalue"]))
        if not r["enough_computed"]:
            L.append("  WARNING: k too small, corrected embedding is short "
                     "({} dims)".format(r["n_corrected_dims"]))
        L.append("  dimensions that are sqrt(deg):  shipped {}/32   corrected {}/{}".format(
            r["shipped_dims_are_degree"], r["corrected_dims_are_degree"],
            r["n_corrected_dims"]))
        L.append("  median R^2 vs sqrt(deg):        shipped {:.4f}   corrected {:.4f}".format(
            r["shipped_r2_median"], r["corrected_r2_median"]))
        L.append("  site rho                        shipped {:.4f}   corrected {:.4f}".format(
            r["site_rho_shipped"], r["site_rho_corrected"]))
        L.append("  kinase rho                      shipped {:.4f}   corrected {:.4f}".format(
            r["kinase_rho_shipped"], r["kinase_rho_corrected"]))
        L.append("")
    txt = "\n".join(L)
    (OUT / "spectral_fix.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
