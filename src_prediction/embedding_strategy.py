"""
Modular embedding strategies for the kinase-substrate LOO pipeline.

To swap node2vec for another graph embedding method, implement
EmbeddingStrategy.fit() and pass an instance to run_leave_one_out().

Current implementations
-----------------------
Node2VecStrategy     : fast DeepWalk (p=q=1) or biased node2vec walks
                       via the custom walk generator in embeddings.py
SpectralEmbeddingStrategy : normalized graph Laplacian eigenvectors (degree-
                       normalizing; addresses node2vec hub bias — Option 1)

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
        from src_prediction.embeddings import fit_node2vec_embeddings

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
    Degree-normalized embedding via the normalized graph Laplacian.

    L_norm = I - D^{-1/2} A D^{-1/2}

    High-degree nodes are penalized by construction: each row/column of A is
    divided by sqrt(degree), so a kinase with 500 graph neighbors does not
    dominate the embedding space merely because of connectivity.

    This directly addresses the hub bias in Node2VecStrategy: ULK1/VRK2/VRK3
    monopolize top-1 predictions in node2vec because they appear in nearly every
    random walk. Spectral embeddings do not share this property.

    Implementation: scipy ARPACK eigsh on L_norm, requesting the `dimensions`
    eigenvectors corresponding to the smallest non-trivial eigenvalues.
    Typical runtime for the liver graph (~16 K nodes, ~2.1 M edges): 60-180 s.

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

        # Request dimensions+1 smallest eigenvectors; drop the trivial λ≈0 one.
        k = min(self.dimensions + 1, n - 1)
        rng = np.random.default_rng(self.random_state)
        eigenvalues, eigenvectors = spla.eigsh(L, k=k, which="SM", v0=rng.standard_normal(n))
        order = np.argsort(eigenvalues)[1 : self.dimensions + 1]
        coords = eigenvectors[:, order]  # (n, dimensions)

        return {str(nodes[i]): coords[i] for i in range(n)}
