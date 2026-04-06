from __future__ import annotations

from src_prediction.config import (
    GENERIC_NODES,
    GENERIC_EDGES,
    LIVER_NODES,
    LIVER_SITE_NODES_OUT,
    CANDIDATE_KINASES_OUT,
    LIVER_POSITIVE_EDGES_OUT,
    PREP_SUMMARY_OUT,
)
from src_prediction.evaluation_set import build_eval_set, summarize_eval_set
from src_prediction.graph_loader import load_graph
from src_prediction.io_utils import write_csv_gz, write_text


def main() -> None:
    print("=== Loading generic graph ===")
    generic = load_graph(GENERIC_NODES, GENERIC_EDGES, label="generic")
    print(f"Generic nodes: {len(generic.nodes):,}")
    print(f"Generic edges: {len(generic.edges):,}\n")

    print("=== Loading liver graph ===")
    liver = load_graph(LIVER_NODES, GENERIC_EDGES, label="liver_nodes_plus_generic_edges_for_prep")
    print(f"Liver nodes: {len(liver.nodes):,}\n")

    print("=== Building evaluation set ===")
    eval_set = build_eval_set(
        generic_nodes=generic.nodes,
        generic_edges=generic.edges,
        liver_nodes=liver.nodes,
    )

    print(f"Liver site nodes: {len(eval_set.liver_site_nodes):,}")
    print(f"Candidate kinases: {len(eval_set.candidate_kinases):,}")
    print(f"Positive kinase-site edges on liver sites: {len(eval_set.positive_edges):,}\n")

    print("=== Writing outputs ===")
    write_csv_gz(eval_set.liver_site_nodes, LIVER_SITE_NODES_OUT)
    write_csv_gz(eval_set.candidate_kinases, CANDIDATE_KINASES_OUT)
    write_csv_gz(eval_set.positive_edges, LIVER_POSITIVE_EDGES_OUT)

    summary = summarize_eval_set(eval_set)
    write_text(summary, PREP_SUMMARY_OUT)

    print(summary)
    print()
    print(f"Wrote: {LIVER_SITE_NODES_OUT}")
    print(f"Wrote: {CANDIDATE_KINASES_OUT}")
    print(f"Wrote: {LIVER_POSITIVE_EDGES_OUT}")
    print(f"Wrote: {PREP_SUMMARY_OUT}")
    print("✅ Done.")


if __name__ == "__main__":
    main()