"""
Run the "predict-all" inference pass (inference_all.run_inference_all) across
all four networks x three embedding strategies = 12 combinations.

Scope (agreed with Dr. Ayati, design discussed 2026-06-23/24):
  Sites:      the 3,228 phosphosites with valid liver FC data (same filter as
              ksea.load_site_fc_values()) -- sites we know are phosphorylated
              in this dataset, not the full 10,146/13,410 liver-site sets.
  Networks:   generic, control, cancer, liver FC (cancer included as a 4th
              heatmap column per Salvador, 2026-06-24 -- Dr. Ayati only asked
              for generic/control/liver FC explicitly, but we already have it).
  Embeddings: node2vec, spectral, and the new ConcatEmbeddingStrategy.

Network graph source: "generic" uses the GENERIC_NODES/GENERIC_EDGES file,
matching the existing LOO scripts (run_multi_kinase_generic*.py) -- NOT the
liver graph restricted to base relations. Stripping correlation edges from
the liver graph fragments it into many small disconnected protein-site
components (liver-only sites whose only edges ARE the correlation edges),
which breaks SpectralEmbeddingStrategy's eigensolver (near-zero eigenvalue
multiplicity blows up ARPACK's unshifted "SM" mode). Consequence: under
"generic", only ~391/3,228 (12%) of the FC-valid sites have any embedding at
all, since the GENERIC graph's site population is essentially the
PSP-annotated sites, a different and much smaller set than the broader
liver-proteomics-detected sites. This is a real, reportable asymmetry (it
underscores how little the generic graph knows about liver-detected
phosphosites), not a bug -- Result 1/2 will show far more missing data for
generic than for the other three networks.

Output: outputs_prediction/inference_all/{network}_{embedding}/
  ranks.csv.gz   -- long format: site_node_id, kinase_node_id, score,
                     rank_in_site, is_true_kinase
  summary.txt    -- coverage + rank distribution stats

Combos that already have a ranks.csv.gz are skipped (safe to re-run after a
partial sweep).
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from src_prediction.core.config import (
    CANDIDATE_KINASES_OUT,
    GENERIC_EDGES,
    GENERIC_NODES,
    LIVER_EDGES,
    LIVER_NODES,
    LIVER_POSITIVE_EDGES_OUT,
    PRED_OUTPUTS_DIR,
)
from src_prediction.embeddings.embedding_strategy import (
    ConcatEmbeddingStrategy,
    Node2VecStrategy,
    SpectralEmbeddingStrategy,
)
from src_prediction.core.graph_loader import load_graph
from src_prediction.engine.inference_all import run_inference_all
from src_prediction.core.io_utils import read_csv_auto, write_csv_gz, write_text
from src_prediction.analysis.ksea import load_site_fc_values
from src_prediction.core.relation_filters import (
    GENERIC_BASE_RELATIONS,
    SITE_CORR_CANCER_RELATIONS,
    SITE_CORR_CTRL_RELATIONS,
    SITE_CORR_FC_RELATIONS,
)

LIVER_EXCEL = Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")
INFERENCE_ALL_DIR = PRED_OUTPUTS_DIR / "inference_all"

NETWORK_CONFIGS = {
    "generic": {
        "nodes": GENERIC_NODES,
        "edges": GENERIC_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS),
    },
    "control": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_CTRL_RELATIONS,
    },
    "cancer": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_CANCER_RELATIONS,
    },
    "liver_fc": {
        "nodes": LIVER_NODES,
        "edges": LIVER_EDGES,
        "allowed_relations": set(GENERIC_BASE_RELATIONS) | SITE_CORR_FC_RELATIONS,
    },
}


def _embedding_strategies() -> dict:
    return {
        "node2vec": Node2VecStrategy(dimensions=32, directed=True, seed=42),
        "spectral": SpectralEmbeddingStrategy(dimensions=32, directed=False, random_state=42),
        "concat": ConcatEmbeddingStrategy(
            node2vec=Node2VecStrategy(dimensions=32, directed=True, seed=42),
            spectral=SpectralEmbeddingStrategy(dimensions=32, directed=False, random_state=42),
            directed=True,
        ),
    }


def _write_summary(out: pd.DataFrame, target_sites: list, network: str, embedding: str) -> str:
    n_sites_total = len(target_sites)
    n_sites_covered = out["site_node_id"].nunique()
    n_true_edges = int(out["is_true_kinase"].sum())
    lines = [
        f"Inference-all summary -- network={network}, embedding={embedding}",
        "=" * 60,
        f"Target sites (FC-valid): {n_sites_total:,}",
        f"Sites with >=1 scored kinase: {n_sites_covered:,} "
        f"({n_sites_covered / n_sites_total:.1%})",
        f"Total (kinase, site) pairs scored: {len(out):,}",
        f"Pairs that are known true KSAs (own-edge masked): {n_true_edges:,}",
        "",
        "Per-kinase average rank (top 20 by lowest mean rank):",
    ]
    avg_rank = (
        out.groupby("kinase_node_id")["rank_in_site"]
        .agg(["mean", "median", "count"])
        .sort_values("mean")
        .head(20)
    )
    for kinase, row in avg_rank.iterrows():
        lines.append(
            f"  {kinase:<20}  mean={row['mean']:.1f}  median={row['median']:.0f}  n={int(row['count'])}"
        )
    return "\n".join(lines)


def main() -> None:
    print("=== Loading FC-valid target sites ===")
    fc_values = load_site_fc_values(LIVER_EXCEL)
    target_sites = sorted(set(fc_values.index.astype(str)))
    print(f"Target sites: {len(target_sites):,}\n")

    print("=== Loading known liver positive edges, restricted to target sites ===")
    all_positive_edges = read_csv_auto(LIVER_POSITIVE_EDGES_OUT)
    target_site_set = set(target_sites)
    positive_edges = all_positive_edges.loc[
        all_positive_edges["site_node_id"].astype(str).isin(target_site_set)
    ].reset_index(drop=True)
    print(f"Restricted positive edges (masking set): {len(positive_edges):,}\n")

    candidate_kinases = read_csv_auto(CANDIDATE_KINASES_OUT)

    embedding_strategies = _embedding_strategies()
    combos = [(net, emb) for net in NETWORK_CONFIGS for emb in embedding_strategies]

    for i, (network, embedding) in enumerate(combos, 1):
        out_dir = INFERENCE_ALL_DIR / f"{network}_{embedding}"
        ranks_path = out_dir / "ranks.csv.gz"
        if ranks_path.exists():
            print(f"[{i}/{len(combos)}] {network} x {embedding} -- already exists, skipping.")
            continue

        print(f"\n{'=' * 70}")
        print(f"[{i}/{len(combos)}] RUNNING: network={network}  embedding={embedding}")
        print("=" * 70)

        cfg = NETWORK_CONFIGS[network]
        graph = load_graph(cfg["nodes"], cfg["edges"], label=network)

        out = run_inference_all(
            graph_edges=graph.edges,
            positive_edges=positive_edges,
            target_sites=target_sites,
            candidate_kinases=candidate_kinases,
            embedding_strategy=embedding_strategies[embedding],
            allowed_relations=cfg["allowed_relations"],
            n_jobs_outer=8,
            max_negatives_per_site=50,
            verbose_every=200,
        )

        write_csv_gz(out, ranks_path)
        summary = _write_summary(out, target_sites, network, embedding)
        write_text(summary, out_dir / "summary.txt")
        print(summary)
        print(f"\nWrote: {ranks_path}")
        print(f"Wrote: {out_dir / 'summary.txt'}")

    print("\nAll combos complete.")


if __name__ == "__main__":
    main()
