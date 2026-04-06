from __future__ import annotations

from typing import Dict, Iterable, Optional, Set, Tuple

import networkx as nx
import numpy as np
import pandas as pd
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

    Returns
    -------
    dict[node_id -> np.ndarray]
    """
    if graph.number_of_nodes() == 0:
        return {}

    n2v = Node2Vec(
        graph,
        dimensions=dimensions,
        walk_length=walk_length,
        num_walks=num_walks,
        workers=workers,
        p=p,
        q=q,
        seed=seed,
    )

    model = n2v.fit(
        window=window,
        min_count=min_count,
        batch_words=batch_words,
    )

    embeddings: Dict[str, np.ndarray] = {}
    for node in graph.nodes():
        key = str(node)
        if key in model.wv:
            embeddings[key] = model.wv[key]

    return embeddings


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