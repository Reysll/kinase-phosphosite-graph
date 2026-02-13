# src/builders/ptmsigdb.py
from __future__ import annotations

import itertools
import re
from typing import List, Optional, Set, Tuple

import pandas as pd

from src.ids import site_id
from src.normalize import normalize_protein, normalize_site
from src.config import MAX_SITE_CLIQUE


# Matches examples like:
#   MAP2K1-S298
#   MAP2K1_S298
#   MAP2K1 Ser298
#   MAP2K1:Ser-298
_SITE_TOKEN_RE = re.compile(
    r"^(?P<gene>[A-Za-z0-9]+)[\s:_\-]*"
    r"(?P<aa>SER|THR|TYR|S|T|Y)[\s:_\-]*"
    r"(?P<pos>\d+)$",
    re.IGNORECASE,
)


def _parse_site_token(token: str) -> Optional[Tuple[str, str]]:
    """
    Parse a token into (gene, site_label), where site_label is like S298.
    Returns None if parsing fails.
    """
    if not token or not isinstance(token, str):
        return None

    tok = token.strip()
    if not tok:
        return None

    m = _SITE_TOKEN_RE.match(tok)
    if not m:
        return None

    gene_raw = m.group("gene")
    aa_raw = m.group("aa")
    pos_raw = m.group("pos")

    gene = normalize_protein(gene_raw)
    site_label = normalize_site(f"{aa_raw}{pos_raw}")
    if site_label is None:
        return None

    return gene, site_label


def build_ptmsigdb_edges(path: str, allowed_site_ids: Set[str]) -> pd.DataFrame:
    """
    Build site-site pathway edges from PTMsigDB.txt.

    Each line is a pathway containing site tokens.
    We connect all SITE nodes that co-occur in the same pathway.

    IMPORTANT:
      allowed_site_ids must contain prefixed IDs like:
        SITE:MAP2K1-S298
    """

    lines_seen = 0
    tokens_seen = 0
    parsed_sites = 0
    kept_sites_in_graph = 0
    rows_with_pairs = 0
    edges_added = 0
    rows_skipped_too_large = 0

    edges: List[Tuple[str, str, str]] = []

    with open(path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            lines_seen += 1
            line = line.strip()
            if not line:
                continue

            parts = re.split(r"[\t\s]+", line)
            if len(parts) < 2:
                continue

            # first column = pathway name (ignored)
            site_tokens = parts[1:]
            tokens_seen += len(site_tokens)

            kept_site_nodes: List[str] = []

            for tok in site_tokens:
                parsed = _parse_site_token(tok)
                if parsed is None:
                    continue

                gene, site_label = parsed
                parsed_sites += 1

                sid = site_id(gene, site_label)  # SITE:GENE-S### 

                if sid in allowed_site_ids:
                    kept_site_nodes.append(sid)
                    kept_sites_in_graph += 1

            kept_site_nodes = sorted(set(kept_site_nodes))
            if len(kept_site_nodes) < 2:
                continue

            if len(kept_site_nodes) > MAX_SITE_CLIQUE:
                rows_skipped_too_large += 1
                continue

            rows_with_pairs += 1

            for a, b in itertools.combinations(kept_site_nodes, 2):
                edges.append((a, b, "site_same_pathway"))
                edges_added += 1

    edges_df = pd.DataFrame(
        edges, columns=["source", "target", "relation"]
    ).drop_duplicates()

    unique_edges = len(edges_df)

    print("PTMsigDB parsing stats:")
    print(f"  lines_seen={lines_seen:,}")
    print(f"  tokens_seen={tokens_seen:,}")
    print(f"  parsed_sites={parsed_sites:,}")
    print(f"  kept_sites_in_graph={kept_sites_in_graph:,}")
    print(f"  rows_with_pairs={rows_with_pairs:,}")
    print(f"  edges_added_raw={edges_added:,}")
    print(f"  edges_unique={unique_edges:,}")
    if rows_skipped_too_large > 0:
        print(
            f"  rows_skipped_too_large={rows_skipped_too_large:,} "
            f"(MAX_SITE_CLIQUE={MAX_SITE_CLIQUE})"
        )

    return edges_df

