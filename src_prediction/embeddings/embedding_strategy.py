"""
Modular embedding strategies for the kinase-substrate LOO pipeline.

To swap node2vec for another graph embedding method, implement
EmbeddingStrategy.fit() and pass an instance to run_leave_one_out().

Current implementations
-----------------------
Node2VecStrategy     : fast DeepWalk (p=q=1) or biased node2vec walks
                       via the custom walk generator in embeddings.py
SpectralEmbeddingStrategy : normalized graph Laplacian eigenvectors. NOTE the
                       `drop_null_space` flag: with the default False the
                       embedding is dominated by sqrt(degree) on a
                       disconnected graph. See the class docstring.
ConcatEmbeddingStrategy   : concatenates Node2Vec + Spectral vectors per node
                       (proposed to Dr. Ayati 2026-06-23, after confirming
                       spectral redistributes the hub bias per-category
                       rather than removing it)

Design note — SpectralEmbeddingStrategy vs. graph contrastive learning (GCL):
  Dr. Ayati suggested GCL with spectral filtering as the degree-normalizing
  option. We implement spectral embedding instead because it achieves the same
  degree-normalization goal with a closed-form solution (no neural network
  training loop, no extra hyperparameters). GCL can replace this via the ABC
  if needed — just implement EmbeddingStrategy.fit().
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict

import networkx as nx
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.sparse.csgraph import connected_components

# Extra eigenvectors requested beyond `dimensions + n_components` when
# drop_null_space=True. ARPACK with which="SM" does not return exactly the k
# smallest eigenvalues on a matrix with a high-multiplicity zero eigenvalue,
# so the buffer absorbs the ones it converges to out of order.
NULL_SPACE_BUFFER = 16


class EmbeddingStrategy(ABC):
    """
    Abstract base for graph embedding strategies.

    Every concrete strategy must implement fit(), which takes a NetworkX
    graph and returns a dict mapping node ID strings to embedding vectors.
    The LOO runner calls fit() exactly once (on the PSP-stripped graph)
    before the trial loop begins.
    """

    # Whether to build a DiGraph (True) or undirected Graph (False).
    # The LOO runner reads this before constructing the embedding graph.
    directed: bool = True

    @abstractmethod
    def fit(self, graph: nx.Graph) -> Dict[str, np.ndarray]:
        """
        Embed every node in *graph*.

        Parameters
        ----------
        graph : nx.Graph or nx.DiGraph

        Returns
        -------
        Dict[str, np.ndarray]
            node_id -> embedding vector (all same length)
        """
        ...


@dataclass
class Node2VecStrategy(EmbeddingStrategy):
    """
    Node2vec / DeepWalk embedding strategy.

    When p=1.0 and q=1.0 (default), uses the fast DeepWalk path that skips
    the O(nodes × avg_degree²) transition-probability precomputation.
    For p≠1 or q≠1, uses biased walks computed on-the-fly (still fast).

    Parameters match fit_node2vec_embeddings() in embeddings.py.
    """

    dimensions: int = 32
    walk_length: int = 10
    num_walks: int = 25
    workers: int = 8
    p: float = 1.0
    q: float = 1.0
    window: int = 5
    min_count: int = 1
    batch_words: int = 4
    seed: int = 42
    directed: bool = True

    def fit(self, graph: nx.Graph) -> Dict[str, np.ndarray]:
        # Local import avoids a circular dependency if embedding_strategy
        # is imported at the top of embeddings.py in the future.
        from src_prediction.embeddings.embeddings import fit_node2vec_embeddings

        return fit_node2vec_embeddings(
            graph=graph,
            dimensions=self.dimensions,
            walk_length=self.walk_length,
            num_walks=self.num_walks,
            workers=self.workers,
            p=self.p,
            q=self.q,
            window=self.window,
            min_count=self.min_count,
            batch_words=self.batch_words,
            seed=self.seed,
        )


@dataclass
class SpectralEmbeddingStrategy(EmbeddingStrategy):
    """
    Embedding from the eigenvectors of the normalized graph Laplacian.

    L_norm = I - D^{-1/2} A D^{-1/2}

    drop_null_space: READ THIS BEFORE USING THE DEFAULT
    ---------------------------------------------------
    The null space of L_norm is spanned by one vector per connected component,
    and on a component that vector is proportional to sqrt(degree) elementwise.
    Its dimension is therefore the number of connected components, not 1.

    With `drop_null_space=False` (the default) fit() requests `dimensions + 1`
    eigenvectors and discards a single one, which is correct only for a
    connected graph. The embedding graphs used in this project are not
    connected: 34 components (generic), 116 (control) and 122 (liver FC) of
    size >= 2. Measured on those graphs, 12 to 17 of the 32 shipped dimensions
    regress on sqrt(degree) at R^2 = 1.0000, and the resulting kinase ranking
    correlates with degree at Spearman rho = -0.9995 (control). So roughly half
    the default embedding is the node's degree rescaled.

    With `drop_null_space=True` fit() requests `dimensions + n_components +
    NULL_SPACE_BUFFER` eigenvectors, keeps the `dimensions` smallest with
    eigenvalue > `null_tol`, and raises RuntimeError if too few clear the
    tolerance. That gives 0 of 32 degree dimensions on the same graphs.

    The default stays False so that results already on disk remain
    reproducible. New work should pass drop_null_space=True.

    What this class does NOT do: it does not correct for hub bias, and earlier
    versions of this docstring claimed it did. The claims that "high-degree
    nodes are penalized by construction" and that spectral embeddings "do not
    share" node2vec's hub bias were wrong in the disconnected case, which is
    every case here. The sqrt(degree) normalization is real for the non-trivial
    eigenvectors, but it does not apply to the null space the default fails to
    discard. Measured against node2vec, the default spectral embedding is the
    MORE degree-dependent of the two (rho -0.9995 against -0.4362). See
    thesis/EMBEDDING_DEFECT.md.

    Implementation: scipy ARPACK eigsh on L_norm with which="SM".
    Typical runtime for the liver graph (~16 K nodes, ~2.1 M edges): 60-180 s
    at the default; requesting the larger k for drop_null_space=True adds to
    that (measure before scheduling a full run).

    Design note — directed vs. undirected:
      directed=False here because the normalized Laplacian is only defined for
      symmetric adjacency matrices. The graph is symmetrized before building A.
      Node2VecStrategy uses directed=True; the LOO runner reads the `directed`
      attribute before constructing the embedding graph, so no manual override
      is needed.
    """

    dimensions: int = 32
    directed: bool = False  # must stay False — spectral embedding requires symmetric adjacency
    random_state: int = 42
    # False reproduces the shipped behaviour, including its degree channel.
    drop_null_space: bool = False
    null_tol: float = 1e-6

    def fit(self, graph: nx.Graph) -> Dict[str, np.ndarray]:
        nodes = list(graph.nodes())
        n = len(nodes)
        if n == 0:
            return {}

        node_to_idx = {nd: i for i, nd in enumerate(nodes)}

        # Build sparse symmetric adjacency (symmetrize directed edges)
        src_list, dst_list = [], []
        for u, v in graph.edges():
            i, j = node_to_idx[u], node_to_idx[v]
            src_list.extend([i, j])
            dst_list.extend([j, i])
        A = sp.csr_matrix(
            (np.ones(len(src_list), dtype=np.float32), (src_list, dst_list)),
            shape=(n, n),
        )
        # Binarize after symmetrization (removes duplicate self-loops)
        A = (A + A.T).sign().astype(np.float32)

        # Normalized Laplacian: L = I - D^{-1/2} A D^{-1/2}
        deg = np.asarray(A.sum(axis=1)).flatten()
        # Isolated nodes (deg=0): set d_inv_sqrt=1 so they get a non-NaN row.
        # Their eigenvector components will cluster near zero and not pollute
        # the informative eigenvectors.
        d_inv_sqrt = sp.diags(np.where(deg > 0, deg**-0.5, 1.0))
        L = sp.eye(n, format="csr") - d_inv_sqrt @ A @ d_inv_sqrt

        if self.drop_null_space:
            # dim(null space) == number of connected components, so ask for
            # enough vectors to clear it and then filter by eigenvalue.
            n_components = connected_components(A, directed=False, return_labels=False)
            k = min(self.dimensions + n_components + NULL_SPACE_BUFFER, n - 1)
        else:
            # Request dimensions+1 smallest eigenvectors; drop the trivial λ≈0 one.
            k = min(self.dimensions + 1, n - 1)
        rng = np.random.default_rng(self.random_state)
        eigenvalues, eigenvectors = spla.eigsh(L, k=k, which="SM", v0=rng.standard_normal(n))

        if self.drop_null_space:
            ranked = np.argsort(eigenvalues)
            informative = [i for i in ranked if eigenvalues[i] > self.null_tol]
            n_dropped = k - len(informative)
            if len(informative) < self.dimensions:
                raise RuntimeError(
                    f"SpectralEmbeddingStrategy: only {len(informative)} of {k} computed "
                    f"eigenvalues exceed null_tol={self.null_tol:g}, need {self.dimensions}. "
                    f"The graph has {n_components} connected components over {n} nodes. "
                    f"Raise NULL_SPACE_BUFFER (currently {NULL_SPACE_BUFFER}), or reduce "
                    f"`dimensions`, or restrict the graph to its giant component."
                )
            order = informative[: self.dimensions]
            print(
                f"[spectral] {n_components} components, {k} eigenvectors computed, "
                f"{n_dropped} null (lambda <= {self.null_tol:g}) dropped, "
                f"keeping {self.dimensions} with lambda in "
                f"[{eigenvalues[order[0]]:.3e}, {eigenvalues[order[-1]]:.3e}]"
            )
        else:
            order = np.argsort(eigenvalues)[1 : self.dimensions + 1]
        coords = eigenvectors[:, order]  # (n, dimensions)

        return {str(nodes[i]): coords[i] for i in range(n)}


@dataclass
class ConcatEmbeddingStrategy(EmbeddingStrategy):
    """
    Concatenates a Node2Vec vector and a Spectral vector per node.

    Motivation (reply to Dr. Ayati, 2026-06-23): spectral embedding fixes the
    GLOBAL hub bias (no single kinase like VRK3 wins everywhere) but does not
    fix the WITHIN-CATEGORY hub bias — one kinase per category (PIK3R1/
    TGFBR2/CDK1) still wins ~100% of the time. The two methods capture
    different structure (node2vec: local random-walk neighbourhoods; spectral:
    global community structure via the Laplacian), so concatenating both
    keeps each signal available to the downstream logistic regression instead
    of forcing one or the other.

    directed vs. undirected:
      The LOO runner builds ONE embedding graph with directedness given by
      this strategy's `directed` attribute, then calls fit(graph) once.
      directed=True here, matching Node2VecStrategy's default, so node2vec
      is fit on exactly the same kind of graph it would see standalone.
      SpectralEmbeddingStrategy.fit() requires symmetric adjacency, so its
      sub-call is given an explicit undirected copy (graph.to_undirected())
      rather than reusing the directed graph node2vec receives — per Dr.
      Ayati's note, each sub-strategy gets its own appropriately-directed
      copy rather than assuming one shared graph works for both.
    """

    node2vec: Node2VecStrategy = field(default_factory=Node2VecStrategy)
    spectral: SpectralEmbeddingStrategy = field(default_factory=SpectralEmbeddingStrategy)
    directed: bool = True

    def fit(self, graph: nx.Graph) -> Dict[str, np.ndarray]:
        n2v_graph = graph if graph.is_directed() else graph.to_directed()
        spectral_graph = graph.to_undirected() if graph.is_directed() else graph

        n2v_emb = self.node2vec.fit(n2v_graph)
        spectral_emb = self.spectral.fit(spectral_graph)

        common_nodes = set(n2v_emb) & set(spectral_emb)
        return {node: np.concatenate([n2v_emb[node], spectral_emb[node]]) for node in common_nodes}
