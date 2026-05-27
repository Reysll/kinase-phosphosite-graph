"""
Modular embedding strategies for the kinase-substrate LOO pipeline.

To swap node2vec for another graph embedding method (e.g. graph contrastive
learning with spectral filtering), implement EmbeddingStrategy.fit() and pass
an instance to run_leave_one_out().

Current implementations
-----------------------
Node2VecStrategy : fast DeepWalk (p=q=1) or biased node2vec walks
                   via the custom walk generator in embeddings.py

Future hooks
------------
GCLStrategy      : graph contrastive learning (not yet implemented)
SpectralStrategy : spectral embedding (not yet implemented)
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Dict

import networkx as nx
import numpy as np


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
