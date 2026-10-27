"""
Test for the `drop_null_space` fix in SpectralEmbeddingStrategy.

WHAT IS BEING TESTED
--------------------
The null space of L_norm = I - D^-1/2 A D^-1/2 is spanned by one vector per
connected component, and restricted to a component that vector is exactly
proportional to sqrt(degree). `fit()` with the shipped default
(drop_null_space=False) requests `dimensions + 1` eigenvectors and discards
ONE, which clears the null space only when the graph is connected. On a graph
with c components it keeps c - 1 null vectors, so that many embedding
dimensions are sqrt(degree) rescaled.

Three checks:

  1. A graph with 3 components. With drop_null_space=False, at least one kept
     dimension regresses on sqrt(degree) at R^2 > 0.99. With True, none do.
  2. A connected graph. Both settings return the same embedding up to sign,
     since there is exactly one null vector and both paths discard it.
  3. A graph with more components than the solver can clear raises
     RuntimeError rather than silently returning degree dimensions.

WHY THE R^2 IS MEASURED WITHIN A COMPONENT
------------------------------------------
ARPACK returns an arbitrary orthonormal basis of the null space, not the
per-component indicator vectors themselves, so a returned null vector is
generally a MIXTURE: a_1 * sqrt(deg) on component 1, a_2 * sqrt(deg) on
component 2, and so on, with different a_i. Regressed against sqrt(degree)
over all nodes at once, such a mixture gives a low R^2 even though every
component-restricted piece is exactly sqrt(degree). Restricting the regression
to one component removes the basis ambiguity and tests what the mathematics
actually guarantees.

The plain all-nodes R^2 is reported alongside, because that is the form used
in the real-graph diagnostic (`thesis/analysis/spectral_fix.py`), where it
reached 1.0000. The two agree there because the real embedding graphs have one
component holding almost every node, so a mixture is dominated by one term. On
a small balanced test graph they do not agree, and only the within-component
form is a valid assertion.

This test does not touch project data or results. Nothing is written.

Usage
-----
    python -m thesis.analysis.test_spectral_nullspace_fix
    python -m pytest thesis/analysis/test_spectral_nullspace_fix.py
"""

from __future__ import annotations

import networkx as nx
import numpy as np

from src_prediction.embeddings.embedding_strategy import SpectralEmbeddingStrategy

DIMS = 4
SEED = 42
R2_THRESHOLD = 0.99


def _r2(y: np.ndarray, x: np.ndarray) -> float:
    """R^2 of y regressed on x with an intercept."""
    X = np.column_stack([x, np.ones(len(y))])
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    resid = y - X @ beta
    ss_tot = float(((y - y.mean()) ** 2).sum())
    if ss_tot == 0.0:
        return float("nan")
    return 1.0 - float((resid**2).sum()) / ss_tot


def _degree_r2_per_dim(graph: nx.Graph, emb: dict) -> tuple[list[float], list[float]]:
    """
    Per dimension, R^2 against sqrt(degree): within the largest component, and
    over all nodes. See the module docstring for why the first is the valid one.
    """
    nodes = list(graph.nodes())
    sqrt_deg = np.sqrt(np.array([graph.degree(nd) for nd in nodes], dtype=float))
    M = np.array([emb[str(nd)] for nd in nodes])

    giant = max(nx.connected_components(graph), key=len)
    mask = np.array([nd in giant for nd in nodes])

    within = [_r2(M[mask, j], sqrt_deg[mask]) for j in range(M.shape[1])]
    overall = [_r2(M[:, j], sqrt_deg) for j in range(M.shape[1])]
    return within, overall


def make_multi_component_graph(sizes=(60, 6, 6), seed=0) -> nx.Graph:
    """Disjoint graphs, each a path plus random chords so degrees vary."""
    rng = np.random.default_rng(seed)
    g = nx.Graph()
    for ci, size in enumerate(sizes):
        nodes = ["c{}_n{}".format(ci, i) for i in range(size)]
        g.add_nodes_from(nodes)
        for a, b in zip(nodes, nodes[1:]):
            g.add_edge(a, b)
        for _ in range(size * 2):
            u, v = rng.choice(size, size=2, replace=False)
            g.add_edge(nodes[u], nodes[v])
    return g


def test_null_space_dimensions_are_degree_by_default() -> None:
    graph = make_multi_component_graph()
    n_comp = nx.number_connected_components(graph)
    assert n_comp == 3, n_comp

    shipped = SpectralEmbeddingStrategy(
        dimensions=DIMS, random_state=SEED, drop_null_space=False).fit(graph)
    fixed = SpectralEmbeddingStrategy(
        dimensions=DIMS, random_state=SEED, drop_null_space=True).fit(graph)

    shipped_within, shipped_all = _degree_r2_per_dim(graph, shipped)
    fixed_within, fixed_all = _degree_r2_per_dim(graph, fixed)

    print("3-component graph: {} nodes, {} components, dimensions={}".format(
        graph.number_of_nodes(), n_comp, DIMS))
    print("  drop_null_space=False")
    print("    R^2 vs sqrt(deg), within giant component : "
          + " ".join("{:.4f}".format(v) for v in shipped_within))
    print("    R^2 vs sqrt(deg), all nodes              : "
          + " ".join("{:.4f}".format(v) for v in shipped_all))
    print("    dims with within-component R^2 > {}: {} of {}".format(
        R2_THRESHOLD, sum(v > R2_THRESHOLD for v in shipped_within), DIMS))
    print("  drop_null_space=True")
    print("    R^2 vs sqrt(deg), within giant component : "
          + " ".join("{:.4f}".format(v) for v in fixed_within))
    print("    R^2 vs sqrt(deg), all nodes              : "
          + " ".join("{:.4f}".format(v) for v in fixed_all))
    print("    dims with within-component R^2 > {}: {} of {}".format(
        R2_THRESHOLD, sum(v > R2_THRESHOLD for v in fixed_within), DIMS))

    n_degree_shipped = sum(v > R2_THRESHOLD for v in shipped_within)
    n_degree_fixed = sum(v > R2_THRESHOLD for v in fixed_within)
    assert n_degree_shipped >= 1, (
        "expected at least one degree dimension with the shipped default, got "
        "{}; within-component R^2 = {}".format(n_degree_shipped, shipped_within))
    assert n_degree_fixed == 0, (
        "expected no degree dimensions with drop_null_space=True, got {}; "
        "within-component R^2 = {}".format(n_degree_fixed, fixed_within))
    # The shipped path keeps c - 1 null vectors out of c.
    assert n_degree_shipped == n_comp - 1, (n_degree_shipped, n_comp)


def test_connected_graph_unchanged_by_the_flag() -> None:
    graph = nx.relabel_nodes(
        nx.barabasi_albert_graph(40, 3, seed=7), lambda i: "n{}".format(i))
    assert nx.number_connected_components(graph) == 1

    shipped = SpectralEmbeddingStrategy(
        dimensions=DIMS, random_state=SEED, drop_null_space=False).fit(graph)
    fixed = SpectralEmbeddingStrategy(
        dimensions=DIMS, random_state=SEED, drop_null_space=True).fit(graph)

    assert set(shipped) == set(fixed)
    nodes = sorted(shipped)
    A = np.array([shipped[nd] for nd in nodes])
    B = np.array([fixed[nd] for nd in nodes])

    cosines = []
    for j in range(DIMS):
        u, v = A[:, j], B[:, j]
        cosines.append(abs(float(u @ v) / (np.linalg.norm(u) * np.linalg.norm(v))))
    print("connected graph: {} nodes, per-dimension |cosine| between settings: ".format(
        graph.number_of_nodes()) + " ".join("{:.6f}".format(c) for c in cosines))

    for j, c in enumerate(cosines):
        assert c > 1.0 - 1e-6, "dimension {} differs between settings, |cosine|={}".format(j, c)


def test_too_many_components_raises() -> None:
    """More null vectors than the solver can clear must fail loudly."""
    graph = make_multi_component_graph(sizes=tuple([2] * 30), seed=1)
    assert nx.number_connected_components(graph) == 30

    try:
        SpectralEmbeddingStrategy(
            dimensions=32, random_state=SEED, drop_null_space=True).fit(graph)
    except RuntimeError as exc:
        print("too-many-components graph raised RuntimeError as intended:")
        print("  {}".format(exc))
        assert "null_tol" in str(exc) and "components" in str(exc)
        return
    raise AssertionError("expected RuntimeError, embedding returned normally")


def main() -> None:
    test_null_space_dimensions_are_degree_by_default()
    print()
    test_connected_graph_unchanged_by_the_flag()
    print()
    test_too_many_components_raises()
    print()
    print("all 3 checks passed")


if __name__ == "__main__":
    main()
