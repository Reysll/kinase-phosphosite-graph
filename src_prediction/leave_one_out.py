from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from src_prediction.embeddings import build_graph_from_edges, fit_node2vec_embeddings
from src_prediction.negative_sampling import sample_negative_pairs_for_site
from src_prediction.pair_features import build_pair_feature_table
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
    held_out_rank = _compute_rank(scored, kinase_node_id=true_kinase)
    best_true_rank = _compute_best_true_rank(scored, true_kinase_ids=true_kinases_for_site)

    top1_pred = scored.iloc[0]["kinase_node_id"] if not scored.empty else None
    top1_score = float(scored.iloc[0]["score"]) if not scored.empty else np.nan

    return {
        "trial_index": trial_index,
        "true_kinase_node_id": true_kinase,
        "site_node_id": site_node,
        "held_out_relation": "phosphorylates",
        "held_out_kinase_rank": held_out_rank,
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


def run_leave_one_out(
    graph_edges: pd.DataFrame,
    positive_edges: pd.DataFrame,
    candidate_kinases: pd.DataFrame,
    node2vec_params: Node2VecParams,
    allowed_relations: Optional[Set[str]] = None,
    max_trials: Optional[int] = None,
    random_state: int = 42,
    verbose_every: int = 10,
    n_jobs_outer: int = 1,
    max_negatives_per_site: int = 50,
    use_threads: bool = True,
) -> pd.DataFrame:
    trials = _prepare_trials(
        positive_edges=positive_edges,
        max_trials=max_trials,
        random_state=random_state,
    )

    candidate_kinase_ids: List[str] = candidate_kinases["node_id"].astype(str).tolist()
    site_to_true_kinases = build_site_to_true_kinases(positive_edges)
    trial_dicts = trials.to_dict(orient="records")

    # --- One-time setup ---

    print("[leave-one-out] Building graph once for all trials...")
    graph = build_graph_from_edges(
        edges_df=graph_edges,
        allowed_relations=allowed_relations,
        directed=node2vec_params.directed,
    )
    print(
        f"[leave-one-out] Graph: {graph.number_of_nodes():,} nodes, "
        f"{graph.number_of_edges():,} edges"
    )

    print("[leave-one-out] Fitting node2vec (once)...")
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
    print(f"[leave-one-out] Embeddings ready: {len(embeddings):,} nodes embedded")

    print("[leave-one-out] Pre-computing training feature matrix (once)...")
    X_train, y_train, train_pair_keys = _build_training_matrix(
        positive_edges=positive_edges,
        embeddings=embeddings,
        candidate_kinase_ids=candidate_kinase_ids,
        site_to_true_kinases=site_to_true_kinases,
        max_negatives_per_site=max_negatives_per_site,
    )
    print(
        f"[leave-one-out] Training matrix: {X_train.shape[0]:,} pairs × {X_train.shape[1]} features"
    )

    unique_test_sites = trials["site_node_id"].astype(str).unique().tolist()
    print(
        f"[leave-one-out] Pre-computing scoring features for {len(unique_test_sites):,} test sites..."
    )
    scoring_cache = _build_scoring_cache(unique_test_sites, candidate_kinase_ids, embeddings)
    print(f"[leave-one-out] Scoring cache ready: {len(scoring_cache):,} sites")

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

    out = run_in_parallel(
        items=tasks,
        worker_fn=_run_single_trial,
        max_workers=n_jobs_outer,
        use_threads=use_threads,
    )

    results = []
    for i, row in enumerate(out, start=1):
        results.append(row)
        if verbose_every and (i % verbose_every == 0 or i == len(tasks)):
            print(f"[leave-one-out] completed {i:,}/{len(tasks):,} trials")

    return pd.DataFrame(results).sort_values("trial_index").reset_index(drop=True)
