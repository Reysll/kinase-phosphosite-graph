from __future__ import annotations

from dataclasses import dataclass

import pandas as pd

from src_prediction.io_utils import read_csv_auto


REQUIRED_NODE_COLS = {
    "node_id",
    "node_type",
    "is_kinase",
    "is_phosphatase",
    "protein_role",
}

REQUIRED_EDGE_COLS = {
    "source",
    "target",
    "relation",
}


@dataclass
class LoadedGraph:
    nodes: pd.DataFrame
    edges: pd.DataFrame


def _validate_nodes(nodes: pd.DataFrame, label: str) -> None:
    missing = REQUIRED_NODE_COLS - set(nodes.columns)
    if missing:
        raise ValueError(f"{label} nodes missing required columns: {sorted(missing)}")


def _validate_edges(edges: pd.DataFrame, label: str) -> None:
    missing = REQUIRED_EDGE_COLS - set(edges.columns)
    if missing:
        raise ValueError(f"{label} edges missing required columns: {sorted(missing)}")


def load_graph(nodes_path, edges_path, label: str) -> LoadedGraph:
    nodes = read_csv_auto(nodes_path)
    edges = read_csv_auto(edges_path)

    _validate_nodes(nodes, label=label)
    _validate_edges(edges, label=label)

    nodes = nodes.copy()
    edges = edges.copy()

    nodes["node_id"] = nodes["node_id"].astype(str)
    nodes["node_type"] = nodes["node_type"].astype(str)
    nodes["protein_role"] = nodes["protein_role"].astype(str)
    nodes["is_kinase"] = nodes["is_kinase"].astype(int)
    nodes["is_phosphatase"] = nodes["is_phosphatase"].astype(int)

    edges["source"] = edges["source"].astype(str)
    edges["target"] = edges["target"].astype(str)
    edges["relation"] = edges["relation"].astype(str)

    return LoadedGraph(nodes=nodes, edges=edges)