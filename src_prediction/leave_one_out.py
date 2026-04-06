from __future__ import annotations

from dataclasses import dataclass
from typing import List, Optional, Set

import numpy as np
import pandas as pd

from src_prediction.embeddings import build_graph_from_edges, fit_node2vec_embeddings
from src_prediction.parallel_utils import run_in_parallel
from src_prediction.scoring import CosineSimilarityScorer


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
    candidate_kinase_ids: List[str]
    node2vec_params: Node2VecParams
    allowed_relations: Optional[Set[str]]


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


def _compute_rank(scored: pd.DataFrame, true_kinase_node_id: str) -> Optional[int]:
    if scored.empty:
        return None

    hit = scored.index[scored["kinase_node_id"] == true_kinase_node_id].tolist()
    if not hit:
        return None

    return int(hit[0] + 1)


def _run_single_fold(task: FoldTask) -> dict:
    fold_row = task.fold_row
    graph_edges = task.graph_edges
    candidate_kinase_ids = task.candidate_kinase_ids
    node2vec_params = task.node2vec_params
    allowed_relations = task.allowed_relations

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

    scorer = CosineSimilarityScorer()

    scored = scorer.score_site_against_kinases(
        site_node_id=site_node,
        candidate_kinases=candidate_kinase_ids,
        embeddings=embeddings,
    )

    rank = _compute_rank(scored, true_kinase_node_id=true_kinase)
    top1_pred = scored.iloc[0]["kinase_node_id"] if not scored.empty else None
    top1_score = float(scored.iloc[0]["score"]) if not scored.empty else np.nan
    num_scored = len(scored)

    return {
        "fold_index": fold_index,
        "true_kinase_node_id": true_kinase,
        "site_node_id": site_node,
        "held_out_relation": relation,
        "rank_of_true_kinase": rank,
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
) -> pd.DataFrame:
    """
    True edge-level leave-one-out with optional parallelization across folds.

    IMPORTANT:
    On Windows, do not combine heavy outer parallelism with inner node2vec
    multiprocessing. Prefer:
      - n_jobs_outer > 1
      - node2vec_params.workers = 1
    """
    folds = _prepare_folds(
        positive_edges=positive_edges,
        max_folds=max_folds,
        random_state=random_state,
    )

    candidate_kinase_ids: List[str] = candidate_kinases["node_id"].astype(str).tolist()
    fold_dicts = folds.to_dict(orient="records")

    tasks = [
        FoldTask(
            fold_row=fold_row,
            graph_edges=graph_edges,
            candidate_kinase_ids=candidate_kinase_ids,
            node2vec_params=node2vec_params,
            allowed_relations=allowed_relations,
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

    results_df = pd.DataFrame(results).sort_values("fold_index").reset_index(drop=True)
    return results_df