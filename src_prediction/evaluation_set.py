from __future__ import annotations

from dataclasses import dataclass

import pandas as pd


@dataclass
class EvalSet:
    liver_site_nodes: pd.DataFrame
    candidate_kinases: pd.DataFrame
    positive_edges: pd.DataFrame


def get_liver_site_nodes(liver_nodes: pd.DataFrame) -> pd.DataFrame:
    """
    Liver-observed phosphosite nodes are site nodes present in the liver network output.
    """
    out = (
        liver_nodes.loc[liver_nodes["node_type"] == "site", ["node_id", "node_type", "protein_role"]]
        .drop_duplicates()
        .reset_index(drop=True)
    )
    return out


def get_candidate_kinases(generic_nodes: pd.DataFrame) -> pd.DataFrame:
    """
    Candidate kinases come from the generic graph metadata.
    We use the generic graph so the candidate set is stable and based on known kinase annotation.
    """
    out = (
        generic_nodes.loc[
            (generic_nodes["node_type"] == "protein") & (generic_nodes["is_kinase"] == 1),
            ["node_id", "node_type", "is_kinase", "is_phosphatase", "protein_role"],
        ]
        .drop_duplicates()
        .sort_values("node_id")
        .reset_index(drop=True)
    )
    return out


def get_liver_positive_kinase_site_edges(
    generic_edges: pd.DataFrame,
    liver_site_nodes: pd.DataFrame,
    candidate_kinases: pd.DataFrame,
) -> pd.DataFrame:
    """
    Known positive kinase-site edges for evaluation:
    - relation must be 'phosphorylates'
    - target must be a liver-observed site
    - source must be one of the candidate kinase proteins
    """
    liver_site_ids = set(liver_site_nodes["node_id"].astype(str))
    kinase_ids = set(candidate_kinases["node_id"].astype(str))

    out = generic_edges.loc[
        (generic_edges["relation"] == "phosphorylates")
        & (generic_edges["target"].isin(liver_site_ids))
        & (generic_edges["source"].isin(kinase_ids)),
        ["source", "target", "relation"],
    ].drop_duplicates()

    out = out.rename(
        columns={
            "source": "kinase_node_id",
            "target": "site_node_id",
        }
    ).reset_index(drop=True)

    return out


def build_eval_set(
    generic_nodes: pd.DataFrame,
    generic_edges: pd.DataFrame,
    liver_nodes: pd.DataFrame,
) -> EvalSet:
    liver_site_nodes = get_liver_site_nodes(liver_nodes)
    candidate_kinases = get_candidate_kinases(generic_nodes)
    positive_edges = get_liver_positive_kinase_site_edges(
        generic_edges=generic_edges,
        liver_site_nodes=liver_site_nodes,
        candidate_kinases=candidate_kinases,
    )

    return EvalSet(
        liver_site_nodes=liver_site_nodes,
        candidate_kinases=candidate_kinases,
        positive_edges=positive_edges,
    )


def summarize_eval_set(eval_set: EvalSet) -> str:
    n_sites = len(eval_set.liver_site_nodes)
    n_kinases = len(eval_set.candidate_kinases)
    n_pos = len(eval_set.positive_edges)
    n_unique_positive_sites = eval_set.positive_edges["site_node_id"].nunique() if n_pos > 0 else 0
    n_unique_positive_kinases = eval_set.positive_edges["kinase_node_id"].nunique() if n_pos > 0 else 0

    lines = [
        "Prediction prep summary",
        "=======================",
        f"Liver site nodes: {n_sites:,}",
        f"Candidate kinase proteins: {n_kinases:,}",
        f"Known positive kinase-site edges on liver sites: {n_pos:,}",
        f"Unique positive liver sites: {n_unique_positive_sites:,}",
        f"Unique positive kinases in evaluation set: {n_unique_positive_kinases:,}",
    ]
    return "\n".join(lines)