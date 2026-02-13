# src/builders/ptmcode2.py
from __future__ import annotations

from collections import Counter
from typing import Dict, List, Optional, Set, Tuple

import pandas as pd

from src.ids import site_id
from src.normalize import normalize_protein, normalize_site


PTMCODE_COLS = [
    "Protein1",
    "Protein2",
    "Species",
    "PTM1",
    "Residue1",
    "rRCS1",
    "Propagated1",
    "PTM2",
    "Residue2",
    "rRCS2",
    "Propagated2",
    "Coevolution_evidence",
    "Manual_evidence",
    "Structure_distance_evidence",
]


def _build_ensp_to_gene_map(gene_map_csv: str) -> Tuple[Dict[str, str], Dict[str, int]]:
    df = pd.read_csv(gene_map_csv)

    expected = {"protein1", "protein2", "gene1", "gene2"}
    if not expected.issubset(set(df.columns)):
        df = df.iloc[:, :4].copy()
        df.columns = ["protein1", "protein2", "gene1", "gene2"]

    counts: Dict[str, Counter] = {}
    rows_seen = 0

    def add_obs(ensp_raw: str, gene_raw: str) -> None:
        ensp = str(ensp_raw).strip()
        gene = str(gene_raw).strip()
        if not ensp or not gene or ensp.lower() == "nan" or gene.lower() == "nan":
            return

        counts.setdefault(ensp, Counter())[gene] += 1

        # store without "9606." prefix too
        if "." in ensp:
            tail = ensp.split(".", 1)[1]
            if tail.startswith("ENSP"):
                counts.setdefault(tail, Counter())[gene] += 1

    for _, row in df.iterrows():
        rows_seen += 1
        add_obs(row["protein1"], row["gene1"])
        add_obs(row["protein2"], row["gene2"])

    ensp_to_gene: Dict[str, str] = {}
    multi_gene = 0

    for ensp, ctr in counts.items():
        gene, _ = ctr.most_common(1)[0]
        ensp_to_gene[ensp] = gene
        if len(ctr) > 1:
            multi_gene += 1

    stats = {
        "map_rows_seen": rows_seen,
        "map_unique_keys": len(ensp_to_gene),
        "map_keys_with_multiple_genes": multi_gene,
    }
    return ensp_to_gene, stats


def _extract_phospho_site_label(residue_raw: str) -> Optional[str]:
    if residue_raw is None:
        return None
    s = str(residue_raw).strip()
    if not s or s.lower() == "nan":
        return None
    # normalize_site returns S298/T468/Y15 or None
    out = normalize_site(s)
    return out


def build_ptmcode2_edges(
    ptmcode_path: str,
    gene_map_csv: str,
    existing_site_ids: Set[str],
    chunksize: int = 250_000,
) -> pd.DataFrame:
    """
    Build SITE-SITE coevolution edges from PTMcode2_associations_between_proteins.txt

    Important: PTMcode2 header is commented out with '##', so we must use skiprows + names.
    """

    ensp_to_gene, map_stats = _build_ensp_to_gene_map(gene_map_csv)

    rows_seen = 0
    rows_human = 0
    rows_phospho_pair = 0
    rows_coevo1 = 0

    rows_dropped_site_parse = 0
    rows_dropped_not_in_graph = 0
    edges_added_raw = 0

    edges: List[Tuple[str, str, str]] = []

    for chunk in pd.read_csv(
        ptmcode_path,
        sep="\t",
        skiprows=4,          # skips the three comment/meta lines at the top
        names=PTMCODE_COLS,  # because the true header line is commented out
        header=None,
        chunksize=chunksize,
        low_memory=False,
    ):
        for _, row in chunk.iterrows():
            rows_seen += 1

            species = str(row["Species"]).strip()
            if species != "Homo sapiens" and species.lower() != "human":
                continue
            rows_human += 1

            ptm1 = str(row["PTM1"]).strip().lower()
            ptm2 = str(row["PTM2"]).strip().lower()
            if ptm1 != "phosphorylation" or ptm2 != "phosphorylation":
                continue
            rows_phospho_pair += 1

            try:
                coevo = int(row["Coevolution_evidence"])
            except Exception:
                continue
            if coevo != 1:
                continue
            rows_coevo1 += 1

            p1_raw = str(row["Protein1"]).strip()
            p2_raw = str(row["Protein2"]).strip()

            # Map Ensembl->gene when possible, else keep as-is (likely already gene)
            g1 = ensp_to_gene.get(p1_raw) or ensp_to_gene.get(p1_raw.split(".", 1)[1] if "." in p1_raw else "")
            g2 = ensp_to_gene.get(p2_raw) or ensp_to_gene.get(p2_raw.split(".", 1)[1] if "." in p2_raw else "")

            if g1 is None:
                g1 = p1_raw
            if g2 is None:
                g2 = p2_raw

            g1 = normalize_protein(g1)
            g2 = normalize_protein(g2)

            site1 = _extract_phospho_site_label(row["Residue1"])
            site2 = _extract_phospho_site_label(row["Residue2"])
            if site1 is None or site2 is None:
                rows_dropped_site_parse += 1
                continue

            sid1 = site_id(g1, site1)
            sid2 = site_id(g2, site2)

            if sid1 not in existing_site_ids or sid2 not in existing_site_ids:
                rows_dropped_not_in_graph += 1
                continue

            a, b = (sid1, sid2) if sid1 <= sid2 else (sid2, sid1)
            edges.append((a, b, "site_coevolution"))
            edges_added_raw += 1

    edges_df = pd.DataFrame(edges, columns=["source", "target", "relation"]).drop_duplicates()

    print("PTMcode2 parsing stats:")
    print(f"  map_rows_seen={map_stats['map_rows_seen']:,}")
    print(f"  map_unique_keys={map_stats['map_unique_keys']:,}")
    print(f"  map_keys_with_multiple_genes={map_stats['map_keys_with_multiple_genes']:,}")
    print(f"  rows_seen={rows_seen:,}")
    print(f"  rows_human={rows_human:,}")
    print(f"  rows_phosphorylation_pairs={rows_phospho_pair:,}")
    print(f"  rows_coevolution_evidence_1={rows_coevo1:,}")
    print(f"  rows_dropped_site_parse={rows_dropped_site_parse:,}")
    print(f"  rows_dropped_not_in_graph={rows_dropped_not_in_graph:,}")
    print(f"  edges_added_raw={edges_added_raw:,}")
    print(f"  edges_unique={len(edges_df):,}")

    return edges_df
