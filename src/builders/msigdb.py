# src/builders/msigdb.py
from __future__ import annotations

import itertools
from typing import List, Set, Tuple

import pandas as pd

from src.ids import protein_id
from src.normalize import normalize_protein
from src.config import MAX_PROTEIN_CLIQUE


def _parse_geneset_line(line: str) -> List[str]:
    """
    Expected format (as you showed):
      GENESET_NAME<TAB>GENE1, GENE2, GENE3, ...

    Returns a list of gene tokens.
    """
    s = line.strip()
    if not s or s.startswith("#"):
        return []

    if "\t" not in s:
        return []

    parts = s.split("\t")
    if len(parts) < 2:
        return []

    gene_list = parts[1].strip()
    if not gene_list:
        return []

    # Split comma-separated gene symbols
    genes = [g.strip() for g in gene_list.split(",") if g.strip()]
    return genes


def build_msigdb_edges(path: str, allowed_protein_ids: Set[str]) -> pd.DataFrame:
    """
    Build protein-protein pathway edges from MsigDB.

    Rule:
      - Each line is a pathway (gene set) with a list of proteins
      - Add edges between all protein pairs in the same pathway
      - Only keep proteins already present in the graph (allowed_protein_ids)
      - Do NOT add new nodes

    allowed_protein_ids must contain prefixed IDs like:
      PROTEIN:MAP2K1
    """
    lines_seen = 0
    tokens_seen = 0
    parsed_proteins = 0
    kept_proteins_in_graph = 0
    rows_with_pairs = 0
    edges_added_raw = 0
    rows_skipped_too_large = 0

    edges: List[Tuple[str, str, str]] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            lines_seen += 1
            gene_tokens = _parse_geneset_line(line)
            if not gene_tokens:
                continue

            tokens_seen += len(gene_tokens)

            kept_nodes: List[str] = []
            for tok in gene_tokens:
                gene = normalize_protein(tok)
                if not gene:
                    continue

                parsed_proteins += 1
                pid = protein_id(gene)

                if pid in allowed_protein_ids:
                    kept_nodes.append(pid)
                    kept_proteins_in_graph += 1

            kept_nodes = sorted(set(kept_nodes))
            if len(kept_nodes) < 2:
                continue

            if len(kept_nodes) > MAX_PROTEIN_CLIQUE:
                rows_skipped_too_large += 1
                continue

            rows_with_pairs += 1
            for a, b in itertools.combinations(kept_nodes, 2):
                edges.append((a, b, "protein_same_pathway"))
                edges_added_raw += 1

    edges_df = pd.DataFrame(edges, columns=["source", "target", "relation"]).drop_duplicates()
    edges_unique = len(edges_df)

    print("MsigDB parsing stats:")
    print(f"  lines_seen={lines_seen:,}")
    print(f"  tokens_seen={tokens_seen:,}")
    print(f"  parsed_proteins={parsed_proteins:,}")
    print(f"  kept_proteins_in_graph={kept_proteins_in_graph:,}")
    print(f"  rows_with_pairs={rows_with_pairs:,}")
    print(f"  edges_added_raw={edges_added_raw:,}")
    print(f"  edges_unique={edges_unique:,}")
    if rows_skipped_too_large > 0:
        print(
            f"  rows_skipped_too_large={rows_skipped_too_large:,} "
            f"(MAX_PROTEIN_CLIQUE={MAX_PROTEIN_CLIQUE})"
        )

    return edges_df
