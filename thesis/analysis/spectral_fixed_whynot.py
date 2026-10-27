"""
Why removing the null space did NOT remove the degree channel.

THE RESULT THAT PROMPTED THIS
-----------------------------
With `drop_null_space=True` the embedding has 0 of 32 dimensions that regress
on sqrt(degree), yet the predicted ranking is still degree-ordered:
Spearman rho went only from -0.9995 to -0.9797 (control) and -0.9973 to -0.9802
(liver FC), and site-invariance got WORSE (eta^2 0.9997 -> 0.9999 / 1.0000,
median per-kinase rank SD 1.05 -> 0.52 and 1.18 -> 0.24). Source:
thesis/analysis/out/spectral_fixed/spectral_fixed_compare.txt.

So the null space was not the binding cause. The candidate explanation is that
degree is in EVERY eigenvector of the normalized Laplacian, not only the null
space. If L_rw = D^-1 (D - A) and L_norm = I - D^-1/2 A D^-1/2, then

    L_norm u = lambda u    <=>    L_rw (D^-1/2 u) = lambda (D^-1/2 u)

so u = D^1/2 v for an eigenvector v of the random-walk Laplacian. Entry i of
ANY L_norm eigenvector therefore carries a factor sqrt(deg_i). The null space
is only the extreme case, where v is constant per component and u is exactly
sqrt(deg). For the non-trivial eigenvectors the sqrt(deg) factor survives as a
per-node SCALE: the embedding vector of node i has norm sqrt(deg_i) * ||v_i||.
Dropping the null space removes the dimensions that are literally degree and
leaves the scale, which is still monotone in degree.

WHAT THIS SCRIPT MEASURES
-------------------------
  1. How many of the 32 kept dimensions regress on sqrt(degree) at R^2 > 0.99,
     under both settings, using the patched class on the real graphs. This
     checks the claim made in the fix commit.

     That R^2 is measured over GIANT-COMPONENT nodes, matching
     `spectral_fix.py`, which is where the 12-to-17-of-32 figure comes from.
     The scope matters: ARPACK returns an arbitrary orthonormal basis of the
     null space, so a returned null vector is a mixture of the per-component
     degree vectors with a different scale on each. Regressed over ALL nodes
     at once such a mixture scores near zero even though every
     component-restricted piece is exactly sqrt(degree), so an all-nodes
     measurement reports 0 of 32 for the shipped setting as well and hides
     the defect entirely. An earlier version of this script made that mistake.
  2. Spearman rho between a kinase's embedding-vector NORM and its degree,
     under both settings. If the explanation above is right this stays near 1
     after the fix, which is the surviving channel.
  3. The same after row-normalizing the embedding to unit length (the
     Ng-Jordan-Weiss convention), which is the candidate remedy.

Whether row normalization changes the RANKING is NOT tested here. That needs
another LOO run. This script only measures the embedding.

Read-only. Writes only into thesis/analysis/out/spectral_fixed/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.spectral_fixed_whynot
"""

from __future__ import annotations

from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd
from scipy import stats

from thesis.analysis.degree_rank import _fit
from thesis.analysis.spectral_nullspace import _embedding_graph

OUT = Path("thesis/analysis/out/spectral_fixed")
DIMS = 32


def _r2_vs_sqrt_deg(M: np.ndarray, sqrt_deg: np.ndarray) -> list[float]:
    return [_fit(sqrt_deg, M[:, j])["pearson_r2"] for j in range(M.shape[1])]


def main() -> None:
    from src_prediction.core.config import CANDIDATE_KINASES_OUT
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy

    OUT.mkdir(parents=True, exist_ok=True)
    cands = read_csv_auto(CANDIDATE_KINASES_OUT)["node_id"].astype(str).tolist()

    rows = []
    for network in ("control", "liver_fc"):
        eg = _embedding_graph(network)
        deg_all = dict(eg.degree())
        giant = max(nx.connected_components(eg), key=len)
        for drop in (False, True):
            emb = SpectralEmbeddingStrategy(
                dimensions=DIMS, directed=False, random_state=42,
                drop_null_space=drop).fit(eg)

            nodes = [nd for nd in eg.nodes() if str(nd) in emb]
            M = np.array([emb[str(nd)] for nd in nodes])
            deg = np.array([deg_all[nd] for nd in nodes], dtype=float)
            sqrt_deg = np.sqrt(deg)

            # giant-component scope for the per-dimension R^2 (see docstring)
            gmask = np.array([nd in giant for nd in nodes])
            r2s = _r2_vs_sqrt_deg(M[gmask], sqrt_deg[gmask])
            r2s_all = _r2_vs_sqrt_deg(M, sqrt_deg)
            norms = np.linalg.norm(M, axis=1)

            # restrict the norm-vs-degree measurement to the scored kinases,
            # which is the population the ranking is computed over
            kin = [i for i, nd in enumerate(nodes) if str(nd) in set(cands)]
            kn, kd = norms[kin], deg[kin]
            rho_kin = stats.spearmanr(kd, kn).statistic
            rho_all = stats.spearmanr(deg, norms).statistic

            # candidate remedy: unit-length rows
            with np.errstate(invalid="ignore", divide="ignore"):
                R = M / norms[:, None]
            R = np.nan_to_num(R)
            r2s_rownorm = _r2_vs_sqrt_deg(R[gmask], sqrt_deg[gmask])

            row = {
                "network": network, "drop_null_space": drop,
                "n_nodes": len(nodes), "n_giant": int(gmask.sum()),
                "n_kinases": len(kin),
                "dims_r2_gt_099": int(sum(v > 0.99 for v in r2s)),
                "median_dim_r2": float(np.median(r2s)),
                "max_dim_r2": float(np.max(r2s)),
                "dims_r2_gt_099_all_nodes_scope": int(sum(v > 0.99 for v in r2s_all)),
                "norm_vs_degree_rho_kinases": float(rho_kin),
                "norm_vs_degree_rho_all_nodes": float(rho_all),
                "norm_vs_sqrt_degree_r2_kinases": _fit(np.sqrt(kd), kn)["pearson_r2"],
                "dims_r2_gt_099_after_rownorm": int(sum(v > 0.99 for v in r2s_rownorm)),
                "median_dim_r2_after_rownorm": float(np.median(r2s_rownorm)),
            }
            rows.append(row)
            print("[whynot] {:<9} drop={!s:<5}: degree dims {:>2}/32, "
                  "norm-vs-degree rho(kinases)={:+.4f}".format(
                      network, drop, row["dims_r2_gt_099"],
                      row["norm_vs_degree_rho_kinases"]))

    df = pd.DataFrame(rows)
    df.to_csv(OUT / "whynot_degree_scale.csv", index=False)

    L = ["WHY THE NULL-SPACE FIX DID NOT REMOVE THE DEGREE CHANNEL",
         "=" * 78, "",
         "Measured on the embedding only. The ranking result it explains is in",
         "spectral_fixed_compare.txt: degree rho moved only -0.9995 -> -0.9797",
         "(control) and -0.9973 -> -0.9802 (liver FC).", ""]
    disp = df.copy()
    for c in ("median_dim_r2", "max_dim_r2", "median_dim_r2_after_rownorm",
              "norm_vs_sqrt_degree_r2_kinases"):
        disp[c] = disp[c].map("{:.4f}".format)
    for c in ("norm_vs_degree_rho_kinases", "norm_vs_degree_rho_all_nodes"):
        disp[c] = disp[c].map("{:+.4f}".format)
    L.append(disp.to_string(index=False))
    L.append("")
    L.append("READING THE COLUMNS")
    L.append("-" * 78)
    L.append("dims_r2_gt_099   how many of the 32 kept dimensions ARE sqrt(degree),")
    L.append("                 over GIANT-COMPONENT nodes. The fix drives this to 0.")
    L.append("..._all_nodes_scope  the same count measured over ALL nodes, which")
    L.append("                 reads 0 under BOTH settings and so hides the defect:")
    L.append("                 ARPACK returns a mixture of the per-component degree")
    L.append("                 vectors, and a mixture does not fit sqrt(degree)")
    L.append("                 globally even though each piece fits it exactly.")
    L.append("norm_vs_degree   Spearman between a kinase's embedding-vector NORM")
    L.append("                 and its graph degree. This is the channel the fix")
    L.append("                 does NOT remove.")
    L.append("after_rownorm    the same dimension count after scaling every row to")
    L.append("                 unit length (Ng-Jordan-Weiss convention).")
    L.append("")
    L.append("THE MECHANISM")
    L.append("-" * 78)
    L.append("L_norm u = lambda u is equivalent to L_rw (D^-1/2 u) = lambda (D^-1/2 u),")
    L.append("so every eigenvector of the normalized Laplacian is u = D^1/2 v for an")
    L.append("eigenvector v of the random-walk Laplacian. Entry i of ANY such")
    L.append("eigenvector carries a factor sqrt(deg_i). The null space is just the")
    L.append("case where v is constant and u is EXACTLY sqrt(deg), which is why")
    L.append("dropping it removes the literal degree dimensions. What survives is")
    L.append("the per-node scale: node i's embedding has norm sqrt(deg_i)*||v_i||.")
    L.append("Pair features [K, S, |K-S|, K*S] are built from those vectors, so a")
    L.append("kinase with large degree still enters the classifier with a larger")
    L.append("vector and keeps its global position.")
    L.append("")
    L.append("WHAT IS NOT TESTED")
    L.append("-" * 78)
    L.append("Whether row normalization changes the held-out ranks. The numbers")
    L.append("above are embedding geometry only. It would need another LOO run per")
    L.append("network, about 1.5 h each at the measured rate.")
    txt = "\n".join(L)
    (OUT / "whynot_degree_scale.txt").write_text(txt, encoding="utf-8")
    print()
    print(txt)


if __name__ == "__main__":
    main()
