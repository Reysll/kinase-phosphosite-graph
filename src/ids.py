# src/ids.py
from __future__ import annotations
from dataclasses import dataclass
from typing import Optional


# -------------------------
# ID constructors (Option B)
# -------------------------

def kinase_id(name: str) -> str:
    return f"KINASE:{name}"


def phosphatase_id(name: str) -> str:
    return f"PHOSPHATASE:{name}"


def protein_id(name: str) -> str:
    return f"PROTEIN:{name}"


def site_id(protein: str, site_label: str) -> str:
    """
    Example:
      site_id("MAP2K1", "S298") -> "SITE:MAP2K1-S298"
    """
    return f"SITE:{protein}-{site_label}"


# -------------------------
# Parsing helpers
# -------------------------

@dataclass(frozen=True)
class ParsedSiteID:
    gene: str
    site_label: str


def parse_site_node_id(node_id: str) -> Optional[ParsedSiteID]:
    """
    Parse SITE-prefixed node IDs.

    Accepts:
      SITE:MAP2K1-S298

    Returns:
      ParsedSiteID(gene="MAP2K1", site_label="S298")
    """
    if not node_id.startswith("SITE:"):
        return None

    body = node_id[len("SITE:"):]
    if "-" not in body:
        return None

    gene, site_label = body.rsplit("-", 1)
    if not gene or not site_label:
        return None

    return ParsedSiteID(gene=gene, site_label=site_label)
