from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Set

import numpy as np
import pandas as pd

from src_prediction.embeddings import build_graph_from_edges, fit_node2vec_embeddings
from src_prediction.model_scoring import train_logistic_pair_model, score_site_with_trained_model
from src_prediction.parallel_utils import run_in_parallel
from src_prediction.truth_utils import build_site_to_true_kinases


@dataclass
class Node2VecParams:
    dimensions: int = 64
    walk_length: int = 30
    num_walks: int = 200
    workers: int = 1
    p: float = 1.0
    q: float = 1.0
    window: int = 10
    min_count: int = 1
    batch_words: int = 4
    seed: int = 42
    directed: bool = True


@dataclass
class FoldTask:
    fold_row: dict
    graph_edges: pd.DataFrame
    all_positive_edges: pd.DataFrame
    candidate_kinase_ids: List[str]
    node2vec_params: Node2VecParams
    allowed_relations: Optional[Set[str]]
    site_to_true_kinases: Dict[str, List[str]]
    max_negatives_per_site: int


def _prepare_folds(
    positive_edges: pd.DataFrame,
    max_folds: Optional[int] = None,
    random_state: int = 42,
) -> pd.DataFrame:
    folds = positive_edges.copy().reset_index(drop=True)

    if max_folds is not None and max_folds < len(folds):
        folds = folds.sample(n=max_folds, random_state=random_state).reset_index(drop=True)

    folds = folds.reset_index(drop=True)
    folds["fold_index"] = np.arange(1, len(folds) + 1)
    return folds


def _compute_rank(scored: pd.DataFrame, kinase_node_id: str) -> Optional[int]:
    if scored.empty:
        return None
    hit = scored.index[scored["kinase_node_id"] == kinase_node_id].tolist()
    if not hit:
        return None
    return int(hit[0] + 1)


def _compute_best_true_rank(scored: pd.DataFrame, true_kinase_ids: List[str]) -> Optional[int]:
    if scored.empty or not true_kinase_ids:
        return None
    true_set = set(map(str, true_kinase_ids))
    hits = scored.index[scored["kinase_node_id"].isin(true_set)].tolist()
    if not hits:
        return None
    return int(min(hits) + 1)


def _run_single_fold(task: FoldTask) -> dict:
    fold_row = task.fold_row
    graph_edges = task.graph_edges
    all_positive_edges = task.all_positive_edges
    candidate_kinase_ids = task.candidate_kinase_ids
    node2vec_params = task.node2vec_params
    allowed_relations = task.allowed_relations
    site_to_true_kinases = task.site_to_true_kinases
    max_negatives_per_site = task.max_negatives_per_site

    fold_index = int(fold_row["fold_index"])
    true_kinase = str(fold_row["kinase_node_id"])
    site_node = str(fold_row["site_node_id"])
    relation = "phosphorylates"

    graph = build_graph_from_edges(
        edges_df=graph_edges,
        drop_edge=(true_kinase, site_node, relation),
        allowed_relations=allowed_relations,
        directed=node2vec_params.directed,
    )

    embeddings = fit_node2vec_embeddings(
        graph=graph,
        dimensions=node2vec_params.dimensions,
        walk_length=node2vec_params.walk_length,
        num_walks=node2vec_params.num_walks,
        workers=node2vec_params.workers,
        p=node2vec_params.p,
        q=node2vec_params.q,
        window=node2vec_params.window,
        min_count=node2vec_params.min_count,
        batch_words=node2vec_params.batch_words,
        seed=node2vec_params.seed,
    )

    train_positive_edges = all_positive_edges.loc[
        ~(
            (all_positive_edges["kinase_node_id"] == true_kinase)
            & (all_positive_edges["site_node_id"] == site_node)
        )
    ].copy()

    trained = train_logistic_pair_model(
        train_positive_edges=train_positive_edges,
        candidate_kinase_ids=candidate_kinase_ids,
        site_to_true_kinases=site_to_true_kinases,
        embeddings=embeddings,
        max_negatives_per_site=max_negatives_per_site,
    )

    scored = score_site_with_trained_model(
        model=trained.model,
        site_node_id=site_node,
        candidate_kinase_ids=candidate_kinase_ids,
        embeddings=embeddings,
    )

    true_kinases_for_site = site_to_true_kinases.get(site_node, [])
    held_out_rank = _compute_rank(scored, kinase_node_id=true_kinase)
    best_true_rank = _compute_best_true_rank(scored, true_kinase_ids=true_kinases_for_site)

    top1_pred = scored.iloc[0]["kinase_node_id"] if not scored.empty else None
    top1_score = float(scored.iloc[0]["score"]) if not scored.empty else np.nan
    num_scored = len(scored)

    return {
        "fold_index": fold_index,
        "true_kinase_node_id": true_kinase,
        "site_node_id": site_node,
        "held_out_relation": relation,
        "held_out_kinase_rank": held_out_rank,
        "best_true_kinase_rank": best_true_rank,
        "n_true_kinases_for_site": len(true_kinases_for_site),
        "all_true_kinases_for_site": ";".join(true_kinases_for_site),
        "top1_predicted_kinase": top1_pred,
        "top1_score": top1_score,
        "num_candidate_kinases_scored": num_scored,
        "site_embedding_present": int(site_node in embeddings),
        "true_kinase_embedding_present": int(true_kinase in embeddings),
    }


def run_leave_one_out(
    graph_edges: pd.DataFrame,
    positive_edges: pd.DataFrame,
    candidate_kinases: pd.DataFrame,
    node2vec_params: Node2VecParams,
    allowed_relations: Optional[Set[str]] = None,
    max_folds: Optional[int] = None,
    random_state: int = 42,
    verbose_every: int = 10,
    n_jobs_outer: int = 1,
    max_negatives_per_site: int = 50,
) -> pd.DataFrame:
    folds = _prepare_folds(
        positive_edges=positive_edges,
        max_folds=max_folds,
        random_state=random_state,
    )

    candidate_kinase_ids: List[str] = candidate_kinases["node_id"].astype(str).tolist()
    site_to_true_kinases = build_site_to_true_kinases(positive_edges)
    fold_dicts = folds.to_dict(orient="records")

    tasks = [
        FoldTask(
            fold_row=fold_row,
            graph_edges=graph_edges,
            all_positive_edges=positive_edges,
            candidate_kinase_ids=candidate_kinase_ids,
            node2vec_params=node2vec_params,
            allowed_relations=allowed_relations,
            site_to_true_kinases=site_to_true_kinases,
            max_negatives_per_site=max_negatives_per_site,
        )
        for fold_row in fold_dicts
    ]

    out = run_in_parallel(
        items=tasks,
        worker_fn=_run_single_fold,
        max_workers=n_jobs_outer,
    )

    results = []
    for i, row in enumerate(out, start=1):
        results.append(row)
        if verbose_every and (i % verbose_every == 0 or i == len(tasks)):
            print(f"[leave-one-out] completed {i:,}/{len(tasks):,} folds")

    return pd.DataFrame(results).sort_values("fold_index").reset_index(drop=True)