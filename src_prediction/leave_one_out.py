from __future__ import annotations

import time
from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from src_prediction.embedding_strategy import EmbeddingStrategy, Node2VecStrategy
from src_prediction.embeddings import build_graph_from_edges
from src_prediction.negative_sampling import sample_negative_pairs_for_site
from src_prediction.pair_features import build_pair_feature_table
from src_prediction.parallel_utils import run_in_parallel
from src_prediction.relation_filters import PSP_KSA_RELATIONS
from src_prediction.truth_utils import build_site_to_true_kinases


# ---------------------------------------------------------------------------
# Backward-compat alias: old callers that imported Node2VecParams from here
# can switch to importing Node2VecStrategy from embedding_strategy instead.
# ---------------------------------------------------------------------------
Node2VecParams = Node2VecStrategy


@dataclass
class TrialTask:
    trial_row: dict
    # Pre-built training matrix — same numpy array reference for every task.
    # With ThreadPoolExecutor (threads share memory) there is zero copying.
    # With fork-based multiprocessing (Linux cluster) it is copy-on-write.
    X_train: np.ndarray  # (N_train, D)
    y_train: np.ndarray  # (N_train,)
    train_pair_keys: np.ndarray  # (N_train,) "kinase_id|site_id" strings
    # Pre-built scoring features for this trial's test site
    site_score_X: Optional[np.ndarray]  # (N_candidates_with_embedding, D)
    site_score_kinase_ids: List[str]
    site_to_true_kinases: Dict[str, List[str]]
    candidate_kinase_ids: List[str]


def _prepare_trials(
    positive_edges: pd.DataFrame,
    max_trials: Optional[int] = None,
    random_state: int = 42,
) -> pd.DataFrame:
    trials = positive_edges.copy().reset_index(drop=True)

    if max_trials is not None and max_trials < len(trials):
        trials = trials.sample(n=max_trials, random_state=random_state).reset_index(drop=True)

    trials = trials.reset_index(drop=True)
    trials["trial_index"] = np.arange(1, len(trials) + 1)
    return trials


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


def _compute_adjusted_held_out_rank(
    scored: pd.DataFrame,
    held_out_kinase: str,
    other_true_kinases: List[str],
) -> Optional[int]:
    """
    Rank of the held-out kinase after dropping all OTHER known true kinases
    from the scored list.

    Motivation: if k2 is another known kinase for this site and is ranked
    above the held-out k1, that is expected correct behaviour, not an error.
    Removing k2 gives a fairer estimate of how hard it is to recover k1.

    Example
    -------
    Scored list: [k7, k4, k2, k1]   (k2 and k1 are both true; k1 is held-out)
    After dropping k2  →  [k7, k4, k1]
    adjusted_held_out_rank = 3   (instead of the naive rank of 4)
    """
    if scored.empty:
        return None
    other_set = set(map(str, other_true_kinases))
    filtered = scored[~scored["kinase_node_id"].isin(other_set)].reset_index(drop=True)
    hit = filtered.index[filtered["kinase_node_id"] == held_out_kinase].tolist()
    if not hit:
        return None
    return int(hit[0] + 1)


def _run_single_trial(task: TrialTask) -> dict:
    trial_row = task.trial_row
    trial_index = int(trial_row["trial_index"])
    true_kinase = str(trial_row["kinase_node_id"])
    site_node = str(trial_row["site_node_id"])

    # Mask out the held-out pair — O(N_train) vectorised numpy op, no Python loop
    held_out_key = f"{true_kinase}|{site_node}"
    mask = task.train_pair_keys != held_out_key
    X = task.X_train[mask]
    y = task.y_train[mask]

    model = LogisticRegression(
        max_iter=1000,
        class_weight="balanced",
        solver="liblinear",
        random_state=42,
    )
    model.fit(X, y)

    # Score using pre-built features for this site
    if task.site_score_X is None or len(task.site_score_X) == 0:
        scored = pd.DataFrame(columns=["kinase_node_id", "site_node_id", "score"])
    else:
        probs = model.predict_proba(task.site_score_X)[:, 1]
        scored = (
            pd.DataFrame(
                {
                    "kinase_node_id": task.site_score_kinase_ids,
                    "site_node_id": site_node,
                    "score": probs,
                }
            )
            .sort_values(["score", "kinase_node_id"], ascending=[False, True])
            .reset_index(drop=True)
        )

    true_kinases_for_site = task.site_to_true_kinases.get(site_node, [])

    # Standard ranks (full scored list)
    held_out_rank = _compute_rank(scored, kinase_node_id=true_kinase)
    best_true_rank = _compute_best_true_rank(scored, true_kinase_ids=true_kinases_for_site)

    # Adjusted held-out rank: drop other true kinases before ranking.
    # This removes the "credit" other co-kinases get from NOT being held out,
    # so the held-out kinase is evaluated only against non-true candidates.
    other_true_kinases = [k for k in true_kinases_for_site if k != true_kinase]
    adjusted_held_out_rank = _compute_adjusted_held_out_rank(
        scored, true_kinase, other_true_kinases
    )

    top1_pred = scored.iloc[0]["kinase_node_id"] if not scored.empty else None
    top1_score = float(scored.iloc[0]["score"]) if not scored.empty else np.nan

    return {
        "trial_index": trial_index,
        "true_kinase_node_id": true_kinase,
        "site_node_id": site_node,
        "held_out_relation": "phosphorylates",
        "held_out_kinase_rank": held_out_rank,
        "adjusted_held_out_rank": adjusted_held_out_rank,
        "best_true_kinase_rank": best_true_rank,
        "n_true_kinases_for_site": len(true_kinases_for_site),
        "all_true_kinases_for_site": ";".join(true_kinases_for_site),
        "top1_predicted_kinase": top1_pred,
        "top1_score": top1_score,
        "num_candidate_kinases_scored": len(scored),
        "site_embedding_present": int(task.site_score_X is not None),
        "true_kinase_embedding_present": int(
            true_kinase in {k for k in task.site_score_kinase_ids}
        ),
    }


def _build_training_matrix(
    positive_edges: pd.DataFrame,
    embeddings: Dict[str, np.ndarray],
    candidate_kinase_ids: List[str],
    site_to_true_kinases: Dict[str, List[str]],
    max_negatives_per_site: int,
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Build the full training feature matrix once.
    Returns (X_train, y_train, train_pair_keys) where train_pair_keys is a
    string array of "kinase_id|site_id" used for per-trial masking.
    """
    neg_tables = []
    for site in sorted(positive_edges["site_node_id"].astype(str).unique()):
        true_set = set(site_to_true_kinases.get(site, []))
        neg_df = sample_negative_pairs_for_site(
            site_node_id=site,
            candidate_kinase_ids=candidate_kinase_ids,
            true_kinases_for_site=true_set,
            max_negatives=max_negatives_per_site,
        )
        neg_tables.append(neg_df)

    neg_all = (
        pd.concat(neg_tables, ignore_index=True)
        if neg_tables
        else pd.DataFrame(columns=["kinase_node_id", "site_node_id", "label"])
    )

    pos_pairs = positive_edges[["kinase_node_id", "site_node_id"]].copy()
    pos_pairs["label"] = 1
    all_pairs = pd.concat([pos_pairs, neg_all], ignore_index=True)

    ft = build_pair_feature_table(all_pairs, embeddings, label_col="label")
    feature_cols = [c for c in ft.columns if c.startswith("f_")]

    X = ft[feature_cols].values
    y = ft["label"].values.astype(int)
    keys = (ft["kinase_node_id"] + "|" + ft["site_node_id"]).values

    return X, y, keys


def _build_scoring_cache(
    sites: List[str],
    candidate_kinase_ids: List[str],
    embeddings: Dict[str, np.ndarray],
) -> Dict[str, Tuple[np.ndarray, List[str]]]:
    """
    Pre-compute scoring feature matrices for every unique test site.
    Returns dict: site_id -> (X_score, kinase_ids_with_embeddings)
    """
    cache: Dict[str, Tuple[np.ndarray, List[str]]] = {}
    for site in sites:
        pairs = pd.DataFrame(
            {
                "kinase_node_id": candidate_kinase_ids,
                "site_node_id": [site] * len(candidate_kinase_ids),
            }
        )
        ft = build_pair_feature_table(pairs, embeddings)
        if not ft.empty:
            feat_cols = [c for c in ft.columns if c.startswith("f_")]
            cache[site] = (ft[feat_cols].values, ft["kinase_node_id"].tolist())
    return cache


def _build_embedding_relations(
    graph_edges: pd.DataFrame,
    allowed_relations: Optional[Set[str]],
) -> Optional[Set[str]]:
    """
    Derive the relation set to use for building the embedding graph.

    PSP KSA edges ('phosphorylates') are excluded because they encode the
    prediction target directly — including them leaks the answer into the
    node2vec embeddings before the LOO trial loop begins.

    If allowed_relations is None (meaning 'use all'), we explicitly enumerate
    all relations present in the edge table and subtract the PSP KSA set.
    """
    if allowed_relations is None:
        all_rels = set(graph_edges["relation"].astype(str).unique())
        return all_rels - PSP_KSA_RELATIONS
    return allowed_relations - PSP_KSA_RELATIONS


def run_leave_one_out(
    graph_edges: pd.DataFrame,
    positive_edges: pd.DataFrame,
    candidate_kinases: pd.DataFrame,
    embedding_strategy: EmbeddingStrategy,
    allowed_relations: Optional[Set[str]] = None,
    max_trials: Optional[int] = None,
    random_state: int = 42,
    verbose_every: int = 10,
    n_jobs_outer: int = 1,
    max_negatives_per_site: int = 50,
    use_threads: bool = True,
    # Legacy keyword kept for scripts that pass node2vec_params= by name
    node2vec_params: Optional[EmbeddingStrategy] = None,
) -> pd.DataFrame:
    """
    Run leave-one-out evaluation.

    Parameters
    ----------
    embedding_strategy : EmbeddingStrategy
        Any concrete EmbeddingStrategy (e.g. Node2VecStrategy).
        PSP 'phosphorylates' edges are removed from the graph before calling
        strategy.fit() to prevent embedding leakage of the prediction target.
    allowed_relations : optional set of str
        Which edge relation types to include in the graph.  None = all.
        PSP KSA edges are additionally excluded at embedding time regardless.
    """
    # Support old callers that pass node2vec_params= kwarg
    if node2vec_params is not None and not isinstance(node2vec_params, EmbeddingStrategy):
        raise TypeError(
            "node2vec_params must be an EmbeddingStrategy instance. "
            "Use Node2VecStrategy from src_prediction.embedding_strategy."
        )
    strategy = node2vec_params if node2vec_params is not None else embedding_strategy

    trials = _prepare_trials(
        positive_edges=positive_edges,
        max_trials=max_trials,
        random_state=random_state,
    )

    candidate_kinase_ids: List[str] = candidate_kinases["node_id"].astype(str).tolist()
    site_to_true_kinases = build_site_to_true_kinases(positive_edges)
    trial_dicts = trials.to_dict(orient="records")

    # --- One-time setup ---
    t_start = time.time()

    # Embedding graph: PSP KSA edges excluded to prevent leakage.
    # The model still TRAINS on those edges as labels (positive_edges),
    # but the embedding only sees structural context (PPI, pathway, coevolution,
    # dephosphorylates, has_site, correlation) not the KSA edges themselves.
    embedding_relations = _build_embedding_relations(graph_edges, allowed_relations)
    print(f"[leave-one-out] Building embedding graph (excluded relations: {PSP_KSA_RELATIONS})...")
    t0 = time.time()
    embedding_graph = build_graph_from_edges(
        edges_df=graph_edges,
        allowed_relations=embedding_relations,
        directed=strategy.directed,
    )
    print(
        f"[leave-one-out] Embedding graph: {embedding_graph.number_of_nodes():,} nodes, "
        f"{embedding_graph.number_of_edges():,} edges  ({time.time() - t0:.1f}s)"
    )

    print(f"[leave-one-out] Fitting embeddings via {type(strategy).__name__} (once)...")
    t0 = time.time()
    embeddings = strategy.fit(embedding_graph)
    print(
        f"[leave-one-out] Embeddings ready: {len(embeddings):,} nodes embedded  ({time.time() - t0:.1f}s)"
    )

    print("[leave-one-out] Pre-computing training feature matrix (once)...")
    t0 = time.time()
    X_train, y_train, train_pair_keys = _build_training_matrix(
        positive_edges=positive_edges,
        embeddings=embeddings,
        candidate_kinase_ids=candidate_kinase_ids,
        site_to_true_kinases=site_to_true_kinases,
        max_negatives_per_site=max_negatives_per_site,
    )
    print(
        f"[leave-one-out] Training matrix: {X_train.shape[0]:,} pairs × {X_train.shape[1]} features"
        f"  ({time.time() - t0:.1f}s)"
    )

    unique_test_sites = trials["site_node_id"].astype(str).unique().tolist()
    print(
        f"[leave-one-out] Pre-computing scoring features for {len(unique_test_sites):,} test sites..."
    )
    t0 = time.time()
    scoring_cache = _build_scoring_cache(unique_test_sites, candidate_kinase_ids, embeddings)
    print(
        f"[leave-one-out] Scoring cache ready: {len(scoring_cache):,} sites  ({time.time() - t0:.1f}s)"
    )

    # --- Build lightweight tasks (all share the same numpy array references) ---
    tasks = []
    for trial_row in trial_dicts:
        site = str(trial_row["site_node_id"])
        score_data = scoring_cache.get(site)
        tasks.append(
            TrialTask(
                trial_row=trial_row,
                X_train=X_train,
                y_train=y_train,
                train_pair_keys=train_pair_keys,
                site_score_X=score_data[0] if score_data else None,
                site_score_kinase_ids=score_data[1] if score_data else [],
                site_to_true_kinases=site_to_true_kinases,
                candidate_kinase_ids=candidate_kinase_ids,
            )
        )

    t_trials_start = time.time()
    n_tasks = len(tasks)

    def _print_progress(done: int, total: int) -> None:
        if verbose_every and (done % verbose_every == 0 or done == total):
            elapsed = time.time() - t_trials_start
            rate = done / elapsed if elapsed > 0 else float("inf")
            eta = (total - done) / rate if rate > 0 else 0.0
            print(
                f"[leave-one-out] {done:,}/{total:,} trials  "
                f"({elapsed:.1f}s elapsed, ETA {eta:.0f}s)",
                flush=True,
            )

    results = run_in_parallel(
        items=tasks,
        worker_fn=_run_single_trial,
        max_workers=n_jobs_outer,
        use_threads=use_threads,
        progress_fn=_print_progress,
    )

    return pd.DataFrame(results).sort_values("trial_index").reset_index(drop=True)
