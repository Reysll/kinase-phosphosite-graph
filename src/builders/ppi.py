# src/builders/ppi.py
from __future__ import annotations

from collections import Counter
from typing import Dict, Set, Tuple, List

import pandas as pd

from src.io_utils import read_table
from src.ids import protein_id


def _build_ensp_to_gene_map(map_path: str) -> Tuple[Dict[str, str], Dict[str, int]]:
    """
    Build a mapping: Ensembl protein id -> gene symbol.

    Expected columns in gene_protein(in).csv:
      protein1, protein2, gene1, gene2

    We will map:
      protein1 -> gene1
      protein2 -> gene2

    If an Ensembl id appears with multiple gene symbols, we keep the most frequent.
    """
    df = read_table(map_path)

    # Try to support both header and no-header cases
    expected_cols = {"protein1", "protein2", "gene1", "gene2"}
    if not expected_cols.issubset(set(df.columns)):
        # fallback to positional
        df = df.iloc[:, :4]
        df.columns = ["protein1", "protein2", "gene1", "gene2"]

    pairs: List[Tuple[str, str]] = []
    for _, row in df.iterrows():
        p1 = str(row["protein1"]).strip()
        p2 = str(row["protein2"]).strip()
        g1 = str(row["gene1"]).strip()
        g2 = str(row["gene2"]).strip()

        if p1 and g1 and p1 != "nan" and g1 != "nan":
            pairs.append((p1, g1))
        if p2 and g2 and p2 != "nan" and g2 != "nan":
            pairs.append((p2, g2))

    # Count most frequent gene per Ensembl id
    counts: Dict[str, Counter] = {}
    for ensp, gene in pairs:
        counts.setdefault(ensp, Counter())[gene] += 1

    ensp_to_gene: Dict[str, str] = {}
    collisions = 0
    for ensp, ctr in counts.items():
        # choose most common gene
        gene, _ = ctr.most_common(1)[0]
        ensp_to_gene[ensp] = gene
        if len(ctr) > 1:
            collisions += 1

    stats = {
        "map_rows_seen": len(df),
        "unique_ensp_mapped": len(ensp_to_gene),
        "ensp_with_multiple_genes": collisions,
    }
    return ensp_to_gene, stats


def build_ppi_edges(
    ppi_path: str,
    gene_map_path: str,
    existing_protein_ids: Set[str],
    min_confidence_score: int = 0,
) -> pd.DataFrame:
    """
    Build protein-protein edges from high_confidence_score(in).csv.

    Expected columns in high_confidence_score(in).csv:
      protein1, protein2, confidence_score

    Returns edges_df with columns:
      source, target, relation, confidence_score
    """
    ensp_to_gene, map_stats = _build_ensp_to_gene_map(gene_map_path)

    ppi_df = read_table(ppi_path)

    expected_cols = {"protein1", "protein2", "confidence_score"}
    if not expected_cols.issubset(set(ppi_df.columns)):
        # fallback to positional
        ppi_df = ppi_df.iloc[:, :3]
        ppi_df.columns = ["protein1", "protein2", "confidence_score"]

    rows_seen = 0
    rows_dropped_score = 0
    rows_dropped_no_map = 0
    rows_dropped_not_in_graph = 0
    edges_added_raw = 0

    edges: List[Tuple[str, str, str, int]] = []

    for _, row in ppi_df.iterrows():
        rows_seen += 1

        p1 = str(row["protein1"]).strip()
        p2 = str(row["protein2"]).strip()

        try:
            score = int(row["confidence_score"])
        except Exception:
            continue

        if score < min_confidence_score:
            rows_dropped_score += 1
            continue

        g1 = ensp_to_gene.get(p1)
        g2 = ensp_to_gene.get(p2)
        if not g1 or not g2:
            rows_dropped_no_map += 1
            continue

        n1 = protein_id(g1)
        n2 = protein_id(g2)

        if n1 not in existing_protein_ids or n2 not in existing_protein_ids:
            rows_dropped_not_in_graph += 1
            continue

        # undirected edge: canonicalize order so duplicates collapse
        a, b = (n1, n2) if n1 <= n2 else (n2, n1)

        edges.append((a, b, "ppi_high_confidence", score))
        edges_added_raw += 1

    edges_df = pd.DataFrame(edges, columns=["source", "target", "relation", "confidence_score"])
    edges_df = edges_df.drop_duplicates(subset=["source", "target", "relation"])

    print("PPI parsing stats:")
    print(f"  lines_seen={rows_seen:,}")
    print(f"  min_confidence_score={min_confidence_score:,}")
    print(f"  mapping: unique_ensp_mapped={map_stats['unique_ensp_mapped']:,}")
    print(f"  mapping: ensp_with_multiple_genes={map_stats['ensp_with_multiple_genes']:,}")
    print(f"  rows_dropped_score={rows_dropped_score:,}")
    print(f"  rows_dropped_no_mapping={rows_dropped_no_map:,}")
    print(f"  rows_dropped_not_in_graph={rows_dropped_not_in_graph:,}")
    print(f"  edges_added_raw={edges_added_raw:,}")
    print(f"  edges_unique={len(edges_df):,}")

    return edges_df
