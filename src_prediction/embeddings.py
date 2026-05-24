from __future__ import annotations

from typing import Dict, List, Optional, Set, Tuple

import networkx as nx
import numpy as np
import pandas as pd
from gensim.models import Word2Vec
from node2vec import Node2Vec


def build_graph_from_edges(
    edges_df: pd.DataFrame,
    drop_edge: Optional[Tuple[str, str, str]] = None,
    allowed_relations: Optional[Set[str]] = None,
    directed: bool = True,
) -> nx.Graph:
    """
    Build a NetworkX graph from edge table.

    Parameters
    ----------
    edges_df : pd.DataFrame
        Must contain source, target, relation.
    drop_edge : optional triple
        Exact (source, target, relation) edge to remove for leave-one-out.
    allowed_relations : optional set[str]
        If provided, only these relation types are kept.
    directed : bool
        Use DiGraph if True, Graph if False.

    Returns
    -------
    nx.Graph or nx.DiGraph
    """
    graph = nx.DiGraph() if directed else nx.Graph()

    work = edges_df.copy()

    if allowed_relations is not None:
        work = work.loc[work["relation"].isin(allowed_relations)].copy()

    if drop_edge is not None:
        src, tgt, rel = drop_edge
        work = work.loc[
            ~(
                (work["source"] == src)
                & (work["target"] == tgt)
                & (work["relation"] == rel)
            )
        ].copy()

    for row in work.itertuples(index=False):
        graph.add_edge(row.source, row.target, relation=row.relation)

    return graph


def _biased_walk_embeddings(
    graph: nx.Graph,
    p: float,
    q: float,
    dimensions: int,
    walk_length: int,
    num_walks: int,
    workers: int,
    window: int,
    min_count: int,
    batch_words: int,
    seed: int,
) -> Dict[str, np.ndarray]:
    """
    Fast biased random walk for p≠1 or q≠1.

    Avoids the node2vec library's O(nodes × avg_degree²) transition-probability
    precomputation by computing per-step weights on-the-fly with integer index
    arrays and numpy (np.isin). Walk-step cost is O(degree) in numpy instead of
    O(degree) in a Python loop, and nothing is precomputed upfront.
    """
    rng = np.random.default_rng(seed)
    nodes = list(graph.nodes())
    node_to_idx = {n: i for i, n in enumerate(nodes)}

    # Adjacency as int32 index arrays — enables fast numpy ops per step
    adj: Dict[int, np.ndarray] = {
        node_to_idx[n]: np.array(
            [node_to_idx[w] for w in graph.neighbors(n)], dtype=np.int32
        )
        for n in nodes
    }

    walks: List[List[str]] = []
    for _ in range(num_walks):
        order = list(range(len(nodes)))
        rng.shuffle(order)
        for start_idx in order:
            walk = [str(nodes[start_idx])]
            cur = start_idx
            prev: Optional[int] = None
            for _ in range(walk_length - 1):
                nbrs = adj[cur]
                if len(nbrs) == 0:
                    break
                if prev is None:
                    nxt = int(nbrs[rng.integers(len(nbrs))])
                else:
                    weights = np.where(
                        nbrs == prev,
                        1.0 / p,
                        np.where(np.isin(nbrs, adj[prev]), 1.0, 1.0 / q),
                    )
                    weights = weights / weights.sum()
                    nxt = int(nbrs[rng.choice(len(nbrs), p=weights)])
                walk.append(str(nodes[nxt]))
                prev = cur
                cur = nxt
            walks.append(walk)

    model = Word2Vec(
        sentences=walks,
        vector_size=dimensions,
        window=window,
        min_count=min_count,
        sg=1,
        workers=workers,
        seed=seed,
        batch_words=batch_words,
    )
    return {str(n): model.wv[str(n)] for n in nodes if str(n) in model.wv}


def _deepwalk_embeddings(
    graph: nx.Graph,
    dimensions: int,
    walk_length: int,
    num_walks: int,
    workers: int,
    window: int,
    min_count: int,
    batch_words: int,
    seed: int,
) -> Dict[str, np.ndarray]:
    """
    Fast DeepWalk for p=q=1: uniform random walks with no transition-probability
    precomputation. Avoids the O(nodes × avg_degree²) bottleneck in the node2vec
    library that dominates runtime on dense graphs.
    """
    rng = np.random.default_rng(seed)
    nodes = list(graph.nodes())

    # Precompute adjacency lists once — O(nodes + edges)
    adj: Dict = {n: list(graph.neighbors(n)) for n in nodes}

    walks: List[List[str]] = []
    for _ in range(num_walks):
        order = nodes.copy()
        rng.shuffle(order)
        for start in order:
            walk = [str(start)]
            cur = start
            for _ in range(walk_length - 1):
                nbrs = adj[cur]
                if not nbrs:
                    break
                cur = nbrs[int(rng.integers(len(nbrs)))]
                walk.append(str(cur))
            walks.append(walk)

    model = Word2Vec(
        sentences=walks,
        vector_size=dimensions,
        window=window,
        min_count=min_count,
        sg=1,
        workers=workers,
        seed=seed,
        batch_words=batch_words,
    )

    return {str(n): model.wv[str(n)] for n in nodes if str(n) in model.wv}


def fit_node2vec_embeddings(
    graph: nx.Graph,
    dimensions: int = 64,
    walk_length: int = 30,
    num_walks: int = 200,
    workers: int = 1,
    p: float = 1.0,
    q: float = 1.0,
    window: int = 10,
    min_count: int = 1,
    batch_words: int = 4,
    seed: int = 42,
) -> Dict[str, np.ndarray]:
    """
    Fit node2vec and return embedding dictionary.
    When p=1 and q=1, uses a fast DeepWalk path that skips node2vec's
    O(nodes × avg_degree²) transition-probability precomputation.
    """
    if graph.number_of_nodes() == 0:
        return {}

    if p == 1.0 and q == 1.0:
        return _deepwalk_embeddings(
            graph=graph,
            dimensions=dimensions,
            walk_length=walk_length,
            num_walks=num_walks,
            workers=workers,
            window=window,
            min_count=min_count,
            batch_words=batch_words,
            seed=seed,
        )

    return _biased_walk_embeddings(
        graph=graph,
        p=p,
        q=q,
        dimensions=dimensions,
        walk_length=walk_length,
        num_walks=num_walks,
        workers=workers,
        window=window,
        min_count=min_count,
        batch_words=batch_words,
        seed=seed,
    )


def embeddings_to_frame(embeddings: Dict[str, np.ndarray]) -> pd.DataFrame:
    """
    Convert embedding dict to a DataFrame for inspection or export.
    """
    if not embeddings:
        return pd.DataFrame(columns=["node_id"])

    rows = []
    for node_id, vec in embeddings.items():
        row = {"node_id": node_id}
        for i, value in enumerate(vec):
            row[f"emb_{i}"] = float(value)
        rows.append(row)

    return pd.DataFrame(rows)