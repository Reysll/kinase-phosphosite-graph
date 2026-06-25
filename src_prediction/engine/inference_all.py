"""
"Predict-all" inference pass — scores every (candidate kinase, target site)
pair without held-out LOO trials, while still avoiding self-leakage for the
known true edges.

Why this exists (agreed with Dr. Ayati, design discussed 2026-06-23):
  LOO only ever evaluates the kinase that was actually held out for a site,
  and only runs on multi-kinase sites (there must be a known edge to remove
  before testing recovery of it). That structurally cannot produce:
    - Result 1: every kinase's rank at every liver-detected site (needs all
      ~420 x ~3,228 scores to build a heatmap/clustergram).
    - Result 2: a full rank distribution per kinase across many sites for a
      KS-test (LOO only gives a kinase's rank at the sites where it happens
      to be a true PSP kinase — too few, too biased a sample).

Two-tier scoring used to avoid leakage while still producing one ranking:
  - Pairs with NO known PSP edge at that site: scored by one globally trained
    model (trained on all known true edges restricted to the target sites,
    nothing held out — there is nothing to leak for a pair that was never a
    training label).
  - Pairs that ARE a known true edge: scored by a model that holds out ONLY
    that one edge (other true edges, including other co-kinases at the same
    site, stay in training). This prevents a kinase's own score from
    trivially reflecting having been a positive training label, without
    penalizing it for co-kinases it has no leakage relationship with.

Reuses the training-matrix and scoring-cache builders from leave_one_out.py
unchanged — only the masking/scoring logic is new.
"""

from __future__ import annotations

import time
from dataclasses import dataclass
from typing import List, Optional, Set

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from src_prediction.embeddings.embedding_strategy import EmbeddingStrategy
from src_prediction.embeddings.embeddings import build_graph_from_edges
from src_prediction.engine.leave_one_out import (
    _build_embedding_relations,
    _build_scoring_cache,
    _build_training_matrix,
)
from src_prediction.core.parallel_utils import run_in_parallel
from src_prediction.data_prep.truth_utils import build_site_to_true_kinases


@dataclass
class MaskedPairTask:
    kinase_node_id: str
    site_node_id: str
    X_train: np.ndarray
    y_train: np.ndarray
    train_pair_keys: np.ndarray
    site_score_X: Optional[np.ndarray]
    site_score_kinase_ids: List[str]


def _fit_model(X: np.ndarray, y: np.ndarray) -> LogisticRegression:
    model = LogisticRegression(
        max_iter=1000,
        class_weight="balanced",
        solver="liblinear",
        random_state=42,
    )
    model.fit(X, y)
    return model


def _score_site(
    model: LogisticRegression,
    site_score_X: Optional[np.ndarray],
    site_score_kinase_ids: List[str],
    site_node_id: str,
) -> pd.DataFrame:
    if site_score_X is None or len(site_score_X) == 0:
        return pd.DataFrame(columns=["kinase_node_id", "site_node_id", "score"])
    probs = model.predict_proba(site_score_X)[:, 1]
    return pd.DataFrame(
        {
            "kinase_node_id": site_score_kinase_ids,
            "site_node_id": site_node_id,
            "score": probs,
        }
    )


def _run_masked_pair(task: MaskedPairTask) -> dict:
    """
    Hold out exactly this one (kinase, site) edge, retrain, and return ONLY
    this kinase's own score at this site. Other kinases' scores from this
    model are intentionally discarded — they get the global model's score
    instead (see module docstring, two-tier scoring).
    """
    held_out_key = f"{task.kinase_node_id}|{task.site_node_id}"
    mask = task.train_pair_keys != held_out_key
    model = _fit_model(task.X_train[mask], task.y_train[mask])

    if task.site_score_X is None or len(task.site_score_X) == 0:
        return {
            "kinase_node_id": task.kinase_node_id,
            "site_node_id": task.site_node_id,
            "score": np.nan,
        }

    try:
        idx = task.site_score_kinase_ids.index(task.kinase_node_id)
    except ValueError:
        return {
            "kinase_node_id": task.kinase_node_id,
            "site_node_id": task.site_node_id,
            "score": np.nan,
        }

    probs = model.predict_proba(task.site_score_X)[:, 1]
    return {
        "kinase_node_id": task.kinase_node_id,
        "site_node_id": task.site_node_id,
        "score": float(probs[idx]),
    }


def run_inference_all(
    graph_edges: pd.DataFrame,
    positive_edges: pd.DataFrame,
    target_sites: List[str],
    candidate_kinases: pd.DataFrame,
    embedding_strategy: EmbeddingStrategy,
    allowed_relations: Optional[Set[str]] = None,
    n_jobs_outer: int = 8,
    max_negatives_per_site: int = 50,
    use_threads: bool = True,
    verbose_every: int = 100,
) -> pd.DataFrame:
    """
    Score every (candidate kinase, target site) pair.

    Parameters
    ----------
    positive_edges : known true kinase-site edges. Used both as training
        labels and as the set of pairs requiring own-edge masking. Should
        already be restricted to target_sites by the caller.
    target_sites : the site set to rank (e.g. the 3,228 FC-valid liver sites).

    Returns
    -------
    Long-format DataFrame, one row per (kinase, site) pair with an embedding
    for both endpoints:
        site_node_id, kinase_node_id, score, rank_in_site, is_true_kinase
    """
    candidate_kinase_ids: List[str] = candidate_kinases["node_id"].astype(str).tolist()
    site_to_true_kinases = build_site_to_true_kinases(positive_edges)

    t0 = time.time()
    embedding_relations = _build_embedding_relations(graph_edges, allowed_relations)
    embedding_graph = build_graph_from_edges(
        edges_df=graph_edges,
        allowed_relations=embedding_relations,
        directed=embedding_strategy.directed,
    )
    print(
        f"[inference-all] Embedding graph: {embedding_graph.number_of_nodes():,} nodes, "
        f"{embedding_graph.number_of_edges():,} edges"
    )

    embeddings = embedding_strategy.fit(embedding_graph)
    print(f"[inference-all] Embeddings ready: {len(embeddings):,} nodes  ({time.time() - t0:.1f}s)")

    t0 = time.time()
    X_train, y_train, train_pair_keys = _build_training_matrix(
        positive_edges=positive_edges,
        embeddings=embeddings,
        candidate_kinase_ids=candidate_kinase_ids,
        site_to_true_kinases=site_to_true_kinases,
        max_negatives_per_site=max_negatives_per_site,
    )
    print(
        f"[inference-all] Training matrix: {X_train.shape[0]:,} x {X_train.shape[1]}"
        f"  ({time.time() - t0:.1f}s)"
    )

    t0 = time.time()
    scoring_cache = _build_scoring_cache(target_sites, candidate_kinase_ids, embeddings)
    print(
        f"[inference-all] Scoring cache: {len(scoring_cache):,}/{len(target_sites):,} sites"
        f"  ({time.time() - t0:.1f}s)"
    )

    print("[inference-all] Fitting global model (no holds)...")
    t0 = time.time()
    global_model = _fit_model(X_train, y_train)
    print(f"[inference-all] Global model fit  ({time.time() - t0:.1f}s)")

    global_rows = []
    for site in target_sites:
        cached = scoring_cache.get(site)
        if cached is None:
            continue
        site_score_X, site_score_kinase_ids = cached
        global_rows.append(_score_site(global_model, site_score_X, site_score_kinase_ids, site))
    global_scores = (
        pd.concat(global_rows, ignore_index=True)
        if global_rows
        else pd.DataFrame(columns=["kinase_node_id", "site_node_id", "score"])
    )
    print(f"[inference-all] Global scores: {len(global_scores):,} pairs")

    # Masked trials: one per known true edge, restricted to target_sites.
    target_site_set = set(target_sites)
    masking_edges = (
        positive_edges.loc[positive_edges["site_node_id"].astype(str).isin(target_site_set)][
            ["kinase_node_id", "site_node_id"]
        ]
        .astype(str)
        .drop_duplicates()
    )

    tasks = []
    for row in masking_edges.itertuples(index=False):
        site = str(row.site_node_id)
        cached = scoring_cache.get(site)
        tasks.append(
            MaskedPairTask(
                kinase_node_id=str(row.kinase_node_id),
                site_node_id=site,
                X_train=X_train,
                y_train=y_train,
                train_pair_keys=train_pair_keys,
                site_score_X=cached[0] if cached else None,
                site_score_kinase_ids=cached[1] if cached else [],
            )
        )

    print(f"[inference-all] Running {len(tasks):,} masked own-edge trials...")
    t0 = time.time()

    def _progress(done: int, total: int) -> None:
        if verbose_every and (done % verbose_every == 0 or done == total):
            print(f"[inference-all] {done:,}/{total:,} masked trials  ({time.time() - t0:.1f}s)")

    masked_results = run_in_parallel(
        items=tasks,
        worker_fn=_run_masked_pair,
        max_workers=n_jobs_outer,
        use_threads=use_threads,
        progress_fn=_progress,
    )
    masked_scores = pd.DataFrame(
        masked_results, columns=["kinase_node_id", "site_node_id", "score"]
    )
    print(f"[inference-all] Masked trials done: {len(masked_scores):,}  ({time.time() - t0:.1f}s)")

    # Assemble final table: global score everywhere, overwritten by the
    # masked own-edge score wherever the pair is a known true edge.
    combined = global_scores.merge(
        masked_scores.rename(columns={"score": "masked_score"}),
        on=["kinase_node_id", "site_node_id"],
        how="left",
    )
    combined["is_true_kinase"] = combined["masked_score"].notna()
    combined["score"] = combined["masked_score"].combine_first(combined["score"])
    combined = combined.drop(columns=["masked_score"])

    combined = combined.sort_values(
        ["site_node_id", "score", "kinase_node_id"], ascending=[True, False, True]
    )
    combined["rank_in_site"] = combined.groupby("site_node_id").cumcount() + 1

    return combined.reset_index(drop=True)
