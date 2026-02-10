# src/builders/ptmsigdb.py
from __future__ import annotations

from itertools import combinations
from typing import List, Set, Tuple

import pandas as pd

from src.normalize import normalize_protein, normalize_site, make_site_id


def _token_to_site_id(token: str) -> str | None:
    """
    Convert a PTMsigDB token into our canonical site id: PROTEIN-S298.

    Accepts tokens like:
      - MAP2K1-S298
      - MAP2K1_S298
      - MAP2K1:S298  (rare, but we sanitize)
    """
    if not token or not isinstance(token, str):
        return None

    t = token.strip().upper()
    if not t:
        return None

    # Remove trailing punctuation that sometimes appears in set files
    t = t.rstrip(",;")

    # Normalize separators like "_" into "-"
    t = t.replace("_", "-")

    if "-" not in t:
        return None

    protein_raw, site_raw = t.rsplit("-", 1)

    protein = normalize_protein(protein_raw)
    site = normalize_site(site_raw)
    if not protein or site is None:
        return None

    return make_site_id(protein, site)


def _parse_line_sites(line: str) -> List[str]:
    """
    PTMsigDB lines may look like:
      PATHWAY_NAME: SITE1 SITE2 SITE3 ...
    or:
      PATHWAY_NAME<TAB>SITE1<TAB>SITE2...

    We return the raw tokens for sites only (not the pathway name).
    """
    s = line.strip()
    if not s or s.startswith("#"):
        return []

    # Split off pathway name if colon exists
    if ":" in s:
        _, rest = s.split(":", 1)
        rest = rest.strip()
        if not rest:
            return []
        # site tokens are whitespace-separated
        return rest.split()

    # Else treat it as tab-delimited or whitespace-delimited:
    # first token is pathway name, remaining are site tokens
    if "\t" in s:
        parts = [p.strip() for p in s.split("\t") if p.strip()]
        return parts[1:] if len(parts) > 1 else []

    parts = s.split()
    return parts[1:] if len(parts) > 1 else []


def build_ptmsigdb_edges(ptmsigdb_path: str, existing_site_ids: Set[str]) -> pd.DataFrame:
    """
    Build site-site edges from PTMsigDB.

    Rule (per your spec):
      - Each line represents one pathway
      - Put edges between ALL sites within the row
      - Ignore sites not already present in our graph
      - Do NOT add new nodes

    Output columns: source, target, relation
    """
    edges: List[Tuple[str, str, str]] = []

    lines_seen = 0
    tokens_seen = 0
    parsed_sites = 0
    kept_sites = 0
    rows_with_pairs = 0

    with open(ptmsigdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            lines_seen += 1
            raw_tokens = _parse_line_sites(line)
            if not raw_tokens:
                continue

            tokens_seen += len(raw_tokens)

            # Convert tokens -> canonical site ids
            site_ids: List[str] = []
            for tok in raw_tokens:
                sid = _token_to_site_id(tok)
                if sid is None:
                    continue
                parsed_sites += 1

                if sid in existing_site_ids:
                    site_ids.append(sid)
                    kept_sites += 1

            # Make unique within row
            site_ids = sorted(set(site_ids))
            if len(site_ids) < 2:
                continue

            rows_with_pairs += 1

            # Add all unordered pairs in the row
            for a, b in combinations(site_ids, 2):
                edges.append((a, b, "site_same_pathway"))

    print(
        "PTMsigDB parsing stats:",
        f"lines_seen={lines_seen:,}",
        f"tokens_seen={tokens_seen:,}",
        f"parsed_sites={parsed_sites:,}",
        f"kept_sites_in_graph={kept_sites:,}",
        f"rows_with_pairs={rows_with_pairs:,}",
        f"edges_added={len(edges):,}",
        sep="\n  ",
    )

    return pd.DataFrame(edges, columns=["source", "target", "relation"])
