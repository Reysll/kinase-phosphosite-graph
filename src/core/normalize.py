# src/normalize.py
from __future__ import annotations
import re

AA_MAP = {
    "SER": "S",
    "THR": "T",
    "TYR": "Y",
    "S": "S",
    "T": "T",
    "Y": "Y",
}

def normalize_protein(name: str) -> str:
    return str(name).strip().upper()

def normalize_site(site_raw: str) -> str | None:
    if site_raw is None:
        return None
    s = str(site_raw).strip().upper()
    if not s:
        return None

    m = re.search(r"(SER|THR|TYR|S|T|Y)[\s\-]*([0-9]+)", s)
    if not m:
        return None

    aa, pos = m.groups()
    aa = AA_MAP.get(aa)
    if aa is None:
        return None

    return f"{aa}{pos}"

def make_site_id(protein: str, site: str) -> str:
    return f"{protein}-{site}"
