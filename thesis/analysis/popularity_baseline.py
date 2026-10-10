"""
Popularity baseline -- rank kinases by substrate count, no embedding.

FAIRNESS PARITY (established by code inspection, see report in out/)
--------------------------------------------------------------------
(a) Training edges. run_multi_kinase_*.py passes
    frozen_trials_multi_kinase.csv.gz as `positive_edges`, and
    run_leave_one_out() uses that SAME frame for _build_training_matrix(),
    build_site_to_true_kinases() and the trial list. The model therefore
    trains on the 6,909 multi-kinase KSA edges only, NOT on all 14,596 PSP
    'phosphorylates' edges. An earlier version of this baseline counted
    substrates over all 14,596, giving it a 2.11x information advantage.
    The PARITY baseline below counts over the 6,909 training edges, with the
    held-out edge removed for its own trial.

(b) Candidate set. Every LOO run reports num_candidate_kinases_scored = 415
    (kinases that have an embedding), not the full 420 candidate list, and the
    scored set is identical across networks. The PARITY baseline ranks those
    same 415.

Both variants are computed and reported so the effect of the fix is visible.

Ordering: substrate count DESC, kinase id ASC -- the tie-break used by the
scoring code. adjusted rank drops other known true kinases for the site first,
matching _compute_adjusted_held_out_rank().

Writes only into thesis/analysis/out/.

Usage
-----
    .venv312/Scripts/python -m thesis.analysis.popularity_baseline
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from thesis.analysis.stat_tests import KEY, RANK_COL, holm, load_all, paired_test, union_unrankable
from thesis.analysis.verify_tests import categories, fmt

PRED = Path("outputs_prediction")
SEED_DIR = Path("thesis/analysis/out/seed_variance")
OUT = Path("thesis/analysis/out")
SEEDS = [42, 7, 13]


def scored_candidate_set() -> list[str]:
    p = PRED / "inference_all" / "liver_fc_spectral" / "ranks.csv.gz"
    return sorted(pd.read_csv(p, usecols=["kinase_node_id"])["kinase_node_id"].unique().tolist())


def build_baseline(keep: set[int], cats: dict[str, str], *,
                   counts_from: str, candidates: list[str]) -> pd.DataFrame:
    """counts_from: 'frozen' (model-parity) or 'psp' (all PSP edges)."""
    from src_prediction.core.config import GENERIC_EDGES, PRED_OUTPUTS_DIR
    from src_prediction.core.io_utils import read_csv_auto
    from src_prediction.data_prep.kinase_categories import compute_psp_substrate_counts
    from src_prediction.data_prep.truth_utils import build_site_to_true_kinases

    trials = read_csv_auto(PRED_OUTPUTS_DIR / "frozen_trials_multi_kinase.csv.gz")
    trials = trials.reset_index(drop=True)
    trials["trial_index"] = np.arange(1, len(trials) + 1)
    site_true = build_site_to_true_kinases(trials)

    if counts_from == "frozen":
        counts = trials["kinase_node_id"].astype(str).value_counts().to_dict()
    else:
        counts = compute_psp_substrate_counts(read_csv_auto(GENERIC_EDGES))

    cand_arr = np.array(candidates)
    idx_of = {k: i for i, k in enumerate(candidates)}
    base_count = np.array([float(counts.get(k, 0)) for k in candidates])
    cat_of = np.array([cats.get(k, "poor") for k in candidates])

    rows = []
    for r in trials.itertuples(index=False):
        ti = int(r.trial_index)
        if ti not in keep:
            continue
        k, s = str(r.kinase_node_id), str(r.site_node_id)
        if k not in idx_of:
            continue
        c = base_count.copy()
        c[idx_of[k]] -= 1.0
        others = {x for x in site_true.get(s, []) if x != k}

        m = np.array([cid not in others for cid in cand_arr])
        ids, cnt = cand_arr[m], c[m]
        srt = ids[np.lexsort((ids, -cnt))]
        pos = int(np.flatnonzero(srt == k)[0]) + 1

        mc = m & (cat_of == cats.get(k, "poor"))
        ids_c, cnt_c = cand_arr[mc], c[mc]
        srt_c = ids_c[np.lexsort((ids_c, -cnt_c))]
        pos_c = int(np.flatnonzero(srt_c == k)[0]) + 1

        rows.append({KEY: ti, "true_kinase_node_id": k, "site_node_id": s,
                     RANK_COL: pos, "rank_in_category": pos_c,
                     "category": cats.get(k, "poor")})
    return pd.DataFrame(rows)


def describe(df: pd.DataFrame, col: str, label: str) -> dict:
    v = df[col].to_numpy(dtype=float)
    return {"group": label, "n": len(v), "median": float(np.median(v)),
            "mean": round(float(v.mean()), 1),
            "Hit@10_pct": round(100.0 * float((v <= 10).mean()), 2),
            "Hit@20_pct": round(100.0 * float((v <= 20).mean()), 2)}


def seed_avg_frame(net: str, keep: set[int]) -> pd.DataFrame | None:
    wide = None
    for s in SEEDS:
        p = (PRED / f"{'generic' if net == 'generic' else 'liver'}_multi_kinase" /
             "results.csv.gz") if s == 42 else \
            SEED_DIR / f"seedvar_{net}_node2vec_seed{s}_results.csv.gz"
        if not p.exists():
            return None
        d = pd.read_csv(p)[[KEY, RANK_COL]]
        d = d[d[KEY].isin(keep)].rename(columns={RANK_COL: f"s{s}"})
        wide = d if wide is None else wide.merge(d, on=KEY)
    cols = [c for c in wide.columns if c != KEY]
    return pd.DataFrame({KEY: wide[KEY], RANK_COL: wide[cols].mean(axis=1)})


def run_tests(base: pd.DataFrame, models: dict[str, pd.DataFrame], keep: set[int],
              *, rank_col: str, subset: pd.Series | None, level: str) -> pd.DataFrame:
    rows = []
    b = base if subset is None else base[subset]
    for label, md in models.items():
        left = (b[[KEY, rank_col, "true_kinase_node_id"]]
                .rename(columns={rank_col: "pop_rank"}))
        right = (md[md[KEY].isin(keep)][[KEY, RANK_COL]]
                 .rename(columns={RANK_COL: "mod_rank"}))
        m = left.merge(right, on=KEY)
        if m.empty:
            continue
        if level == "kinase":
            g = m.groupby("true_kinase_node_id")[["pop_rank", "mod_rank"]].median()
            a, bb = g["pop_rank"].to_numpy(float), g["mod_rank"].to_numpy(float)
        else:
            a, bb = m["pop_rank"].to_numpy(float), m["mod_rank"].to_numpy(float)
        if len(a) < 10:
            continue
        r = paired_test(a, bb, "Popularity", label)
        d = a - bb
        r["B_wins"] = int((d > 0).sum()); r["A_wins"] = int((d < 0).sum())
        r["ties"] = int((d == 0).sum())
        rows.append(r)
    if not rows:
        return pd.DataFrame()
    t = pd.DataFrame(rows)
    t["p_holm"] = holm(t["p_value"].tolist())
    t["significant_holm_0.05"] = t["p_holm"] < 0.05
    return t[["comparison", "n", "median_A", "median_B", "median_paired_diff",
              "rank_biserial", "B_wins", "A_wins", "ties", "p_value", "p_holm",
              "significant_holm_0.05", "favours"]]


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    frames, _ = load_all()
    keep = set(frames[("Generic", "node2vec")][KEY].astype(int)) - union_unrankable(frames)
    cats = categories()
    scored = scored_candidate_set()

    L: list[str] = ["POPULARITY BASELINE -- PARITY-CORRECTED", "=" * 78, ""]
    add = L.append
    add("(a) Model trains on the 6,909 frozen multi-kinase edges, NOT all 14,596")
    add("    PSP edges. Code path: run_multi_kinase_*.py passes")
    add("    frozen_trials_multi_kinase.csv.gz as positive_edges ->")
    add("    run_leave_one_out -> _build_training_matrix(positive_edges=...).")
    add("(b) Models score 415 candidates (num_candidate_kinases_scored = 415 in all")
    add("    four runs); the scored set is identical across networks.")
    add("")

    par = build_baseline(keep, cats, counts_from="frozen", candidates=scored)
    old = build_baseline(keep, cats, counts_from="psp",
                         candidates=sorted(set(scored) | set(cats)))

    add("EFFECT OF THE PARITY FIX (global adjusted rank, all trials)")
    add("-" * 78)
    cmpt = pd.DataFrame([
        {**describe(old, RANK_COL, "OLD: all 14,596 PSP counts, 420 candidates")},
        {**describe(par, RANK_COL, "PARITY: 6,909 training counts, 415 candidates")},
    ])
    add(cmpt.to_string(index=False))
    add("")

    for col, title in ((RANK_COL, "GLOBAL"), ("rank_in_category", "WITHIN-CATEGORY")):
        add(f"{title} ADJUSTED RANK -- parity baseline")
        add("-" * 78)
        rows = [describe(par, col, "ALL")]
        for c in ("rich", "average", "poor"):
            sub = par[par["category"] == c]
            if len(sub):
                rows.append(describe(sub, col, c))
        tt = pd.DataFrame(rows)
        add(tt.to_string(index=False))
        tt.to_csv(OUT / f"popularity_parity_{title.lower()}.csv", index=False)
        add("")

    models = {
        "Generic spectral": frames[("Generic", "spectral")],
        "Liver FC spectral": frames[("Liver FC", "spectral")],
    }
    for net, key in (("generic", "Generic"), ("liver_fc", "Liver FC")):
        sa = seed_avg_frame(net, keep)
        if sa is not None:
            models[f"{key} n2v seed-avg"] = sa

    add("PAIRED COMPARISONS -- ALL TRIALS (median_A = popularity, B = model)")
    add("-" * 78)
    t_all = run_tests(par, models, keep, rank_col=RANK_COL, subset=None, level="trial")
    add(fmt(t_all))
    t_all.to_csv(OUT / "popularity_parity_tests_all.csv", index=False)
    add("")

    add("STRATIFIED COMPARISONS -- does any model beat popularity where")
    add("annotation is sparse?")
    add("=" * 78)
    out_frames = []
    for stratum in ("poor", "average"):
        sub = par["category"] == stratum
        for col, cname in ((RANK_COL, "global"), ("rank_in_category", "within-cat")):
            for level in ("trial", "kinase"):
                t = run_tests(par, models, keep, rank_col=col, subset=sub, level=level)
                if t.empty:
                    continue
                add(f"--- {stratum.upper()} | {cname} rank | {level} level ---")
                add(fmt(t))
                add("")
                out_frames.append(t.assign(stratum=stratum, rank_type=cname, level=level))
    if out_frames:
        pd.concat(out_frames, ignore_index=True).to_csv(
            OUT / "popularity_parity_stratified.csv", index=False)

    par.to_csv(OUT / "popularity_parity_trials.csv.gz", index=False, compression="gzip")
    text = "\n".join(L)
    (OUT / "popularity_baseline.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
