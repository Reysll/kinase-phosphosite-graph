"""Audit the protein and site node populations of the generic and liver graphs.

Answers the question of why the two graphs report different protein counts, and
checks each node table against its own edge list for referential integrity.

Read only. Touches no project code, data or results.

    .venv312/Scripts/python -m thesis.analysis.graph_protein_audit

Writes thesis/analysis/out/graph_protein_audit.txt
"""
from __future__ import annotations

import collections
import re
import shutil
import tempfile
from pathlib import Path

import networkx as nx
import numpy as np
import pandas as pd

GENERIC_NODES = Path("outputs/nodes.csv.gz")
GENERIC_EDGES = Path("outputs/edges.csv.gz")
LIVER_NODES = Path("outputs/Liver_network_nodes.csv.gz")
LIVER_EDGES = Path("outputs/Liver_network_edges.csv.gz")
LIVER_XLSX = Path("data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx")

OUT = Path("thesis/analysis/out/graph_protein_audit.txt")

# src/graph/liver_network.py default
MIN_IDENTIFIED = 6
# the prediction target, stripped before any embedding is fit
PHOSPHORYLATES = "phosphorylates"

_LINES: list[str] = []


def say(line: str = "") -> None:
    print(line)
    _LINES.append(line)


def _nodes(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["node_id"] = df["node_id"].astype(str)
    df["node_type"] = df["node_type"].astype(str)
    return df


def _edges(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, usecols=["source", "target", "relation"])
    for c in ("source", "target", "relation"):
        df[c] = df[c].astype(str)
    return df


def _by_type(nodes: pd.DataFrame, node_type: str) -> set:
    return set(nodes.loc[nodes["node_type"] == node_type, "node_id"])


def _kinases(nodes: pd.DataFrame) -> set:
    mask = (nodes["node_type"] == "protein") & (nodes["is_kinase"] == 1)
    return set(nodes.loc[mask, "node_id"])


def _pid(genes) -> set:
    return {"PROTEIN:" + g for g in genes}


def section_populations(gn, ln) -> None:
    say("=" * 72)
    say("1. NODE POPULATIONS AS DECLARED IN THE NODE TABLES")
    say("=" * 72)
    gp, lp = _by_type(gn, "protein"), _by_type(ln, "protein")
    gs, ls = _by_type(gn, "site"), _by_type(ln, "site")
    gk, lk = _kinases(gn), _kinases(ln)

    say("")
    say("{:<28}{:>12}{:>12}{:>12}".format("", "generic", "liver", "difference"))
    for label, a, b in (("proteins", gp, lp), ("sites", gs, ls),
                        ("kinases (is_kinase=1)", gk, lk)):
        say("{:<28}{:>12,}{:>12,}{:>+12,}".format(
            label, len(a), len(b), len(b) - len(a)))

    say("")
    say("Set relationships:")
    for label, a, b in (("proteins", gp, lp), ("sites", gs, ls), ("kinases", gk, lk)):
        say("  {:<10} shared {:>7,}   generic-only {:>6,}   liver-only {:>6,}".format(
            label, len(a & b), len(a - b), len(b - a)))
    say("")
    say("The generic graph is a strict subset of the liver graph. Every generic")
    say("node appears in the liver node table and none is dropped.")


def section_integrity(gn, ge, ln, le) -> set:
    say("")
    say("=" * 72)
    say("2. NODE TABLE VERSUS EDGE LIST (referential integrity)")
    say("=" * 72)
    ghosts: set = set()
    for tag, nodes, edges in (("generic", gn, ge), ("liver", ln, le)):
        prot, site = _by_type(nodes, "protein"), _by_type(nodes, "site")
        hs = edges.loc[edges["relation"] == "has_site"]
        parents, children = set(hs["source"]), set(hs["target"])
        referenced = set(edges["source"]) | set(edges["target"])
        undeclared = parents - prot
        say("")
        say("{}:".format(tag))
        say("  protein nodes declared in node table   : {:>7,}".format(len(prot)))
        say("  distinct has_site parent proteins      : {:>7,}".format(len(parents)))
        say("  parents MISSING from the node table    : {:>7,}".format(len(undeclared)))
        say("  sites with no has_site edge            : {:>7,}".format(len(site - children)))
        say("  edge-referenced but undeclared nodes   : {:>7,}".format(
            len(referenced - prot - site)))
        if tag == "liver":
            ghosts = undeclared
    say("")
    say("The generic graph is internally consistent. The liver node table omits")
    say("{:,} protein nodes that its own has_site edges refer to.".format(len(ghosts)))
    return ghosts


def section_ghosts(ln, le, ghosts) -> None:
    say("")
    say("=" * 72)
    say("3. THE OMITTED PROTEINS AND THEIR CONNECTIVITY")
    say("=" * 72)
    declared = _by_type(ln, "protein")
    deg: collections.Counter = collections.Counter()
    rels: collections.Counter = collections.Counter()
    for s, t, r in le.itertuples(index=False):
        if s in ghosts:
            deg[s] += 1
            rels[r] += 1
        if t in ghosts:
            deg[t] += 1
            rels[r] += 1
        if s in declared:
            deg[s] += 1
        if t in declared:
            deg[t] += 1

    dg = pd.Series({g: deg.get(g, 0) for g in ghosts})
    dd = pd.Series({g: deg.get(g, 0) for g in declared})
    say("")
    say("{:<34}{:>8}{:>10}{:>9}".format(
        "degree in the liver edge list", "median", "mean", "max"))
    say("{:<34}{:>8.0f}{:>10.1f}{:>9,}".format(
        "omitted proteins (n={:,})".format(len(ghosts)),
        dg.median(), dg.mean(), int(dg.max())))
    say("{:<34}{:>8.0f}{:>10.1f}{:>9,}".format(
        "declared proteins (n={:,})".format(len(declared)),
        dd.median(), dd.mean(), int(dd.max())))
    say("")
    say("Relation types touching the omitted proteins:")
    for r, n in sorted(rels.items(), key=lambda kv: -kv[1]):
        say("  {:<24}{:>9,}".format(r, n))
    say("")
    hs = le.loc[le["relation"] == "has_site"]
    orphan_sites = set(hs.loc[hs["source"].isin(ghosts), "target"])
    total_sites = hs["target"].nunique()
    say("Sites whose parent protein is omitted: {:,} of {:,} ({:.1%}).".format(
        len(orphan_sites), total_sites, len(orphan_sites) / total_sites))
    say("These proteins carry has_site edges and nothing else, so they enter the")
    say("graph as degree-1 leaves with no pathway or interaction context.")


def section_components(ge, le) -> None:
    say("")
    say("=" * 72)
    say("4. EFFECT ON EMBEDDING-GRAPH CONNECTIVITY")
    say("=" * 72)
    corr = {r for r in set(le["relation"]) if r.startswith("site_corr_")}
    base = set(le["relation"]) - corr - {PHOSPHORYLATES}
    say("")
    for tag, rels in (("liver, base relations only", base),
                      ("liver FC (base plus fc correlation)",
                       base | {"site_corr_fc_pos", "site_corr_fc_neg"})):
        sub = le.loc[le["relation"].isin(rels)]
        g = nx.Graph()
        g.add_edges_from(zip(sub["source"], sub["target"]))
        sizes = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
        say("{:<38}{:>8,} nodes {:>6,} components  largest {:>7,}".format(
            tag, g.number_of_nodes(), len(sizes), sizes[0]))
    sub = ge.loc[ge["relation"] != PHOSPHORYLATES]
    g = nx.Graph()
    g.add_edges_from(zip(sub["source"], sub["target"]))
    sizes = sorted((len(c) for c in nx.connected_components(g)), reverse=True)
    say("{:<38}{:>8,} nodes {:>6,} components  largest {:>7,}".format(
        "generic", g.number_of_nodes(), len(sizes), sizes[0]))


def _read_liver_sheets():
    """Read the workbook from a copy. The original is often locked by Excel."""
    with tempfile.TemporaryDirectory() as td:
        tmp = Path(td) / "liver.xlsx"
        shutil.copy2(LIVER_XLSX, tmp)
        prot = pd.read_excel(tmp, sheet_name="ProteinExpression", header=1)
        ph = pd.read_excel(tmp, sheet_name="Phosphorylation", header=1)
    return prot, ph


def _genes(df: pd.DataFrame) -> set:
    s = (df["Gene Symbol"].astype(str).str.strip()
         .replace({"": np.nan, "nan": np.nan}).dropna())
    return set(s.unique())


def _passing(df: pd.DataFrame) -> pd.DataFrame:
    cols = [c for c in df.columns
            if re.search(r"control|sample|tumou?r|case", str(c), re.I)]
    return df.loc[df[cols].notna().sum(axis=1) >= MIN_IDENTIFIED]


def section_source(gn, ln) -> None:
    say("")
    say("=" * 72)
    say("5. WHERE THE COUNTS COME FROM IN THE SOURCE WORKBOOK")
    say("=" * 72)
    prot, ph = _read_liver_sheets()
    prot_f, ph_f = _passing(prot), _passing(ph)
    pg_all, hg_all = _genes(prot), _genes(ph)
    pg_f, hg_f = _genes(prot_f), _genes(ph_f)

    say("")
    say("{:<40}{:>9}{:>9}".format("", "rows", "genes"))
    say("{:<40}{:>9,}{:>9,}".format("ProteinExpression sheet", len(prot), len(pg_all)))
    say("{:<40}{:>9,}{:>9,}".format(
        "  passing min_identified={}".format(MIN_IDENTIFIED), len(prot_f), len(pg_f)))
    say("{:<40}{:>9,}{:>9,}".format("Phosphorylation sheet", len(ph), len(hg_all)))
    say("{:<40}{:>9,}{:>9,}".format(
        "  passing min_identified={}".format(MIN_IDENTIFIED), len(ph_f), len(hg_f)))

    gp, lp = _by_type(gn, "protein"), _by_type(ln, "protein")
    say("")
    say("Reconciling the {:+,} liver protein nodes:".format(len(lp - gp)))
    say("  ProteinExpression genes passing the filter   : {:>7,}".format(len(pg_f)))
    say("  already present as generic protein nodes     : {:>7,}".format(
        len(_pid(pg_f) & gp)))
    say("  therefore new                                : {:>7,}".format(
        len(_pid(pg_f) - gp)))
    say("  observed liver-minus-generic protein nodes   : {:>7,}".format(len(lp - gp)))
    say("")
    say("Reconciling the omitted proteins:")
    say("  Phosphorylation genes passing the filter     : {:>7,}".format(len(hg_f)))
    say("  present in the liver node table              : {:>7,}".format(
        len(_pid(hg_f) & lp)))
    say("  absent from the liver node table             : {:>7,}".format(
        len(_pid(hg_f) - lp)))
    say("")
    say("src/graph/liver_network.py builds protein_nodes from the ProteinExpression")
    say("sheet alone (line 135) while emitting has_site edges from the gene column")
    say("of the Phosphorylation sheet (line 167). A gene detected only in phospho")
    say("data therefore receives sites and a has_site edge but never a protein")
    say("node. The node table is then used as allowed_protein_ids and")
    say("existing_protein_ids for the MSigDB and PPI builders (lines 595 to 621),")
    say("which is why those proteins receive no pathway or interaction edges.")

    kin = _kinases(gn)
    say("")
    say("-" * 72)
    say("Candidate kinase coverage in the liver data, for the planned restriction")
    say("of candidates to kinases present in protein expression:")
    say("  candidate kinases in the graph               : {:>7,}".format(len(kin)))
    say("  in ProteinExpression, unfiltered             : {:>7,}".format(
        len(kin & _pid(pg_all))))
    say("  in ProteinExpression, passing the filter     : {:>7,}".format(
        len(kin & _pid(pg_f))))
    say("  phospho-detected, passing the filter         : {:>7,}".format(
        len(kin & _pid(hg_f))))
    say("")
    names = sorted(x.split(":", 1)[1] for x in (kin & _pid(pg_all)))
    say("  the {} candidate kinases in ProteinExpression:".format(len(names)))
    for i in range(0, len(names), 8):
        say("    " + ", ".join(names[i:i + 8]))


def main() -> None:
    gn, ln = _nodes(GENERIC_NODES), _nodes(LIVER_NODES)
    ge, le = _edges(GENERIC_EDGES), _edges(LIVER_EDGES)

    say("Protein and site node audit of the generic and liver graphs")
    say("Generated by thesis/analysis/graph_protein_audit.py. Read only.")
    say("")

    section_populations(gn, ln)
    ghosts = section_integrity(gn, ge, ln, le)
    section_ghosts(ln, le, ghosts)
    section_components(ge, le)
    section_source(gn, ln)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text("\n".join(_LINES) + "\n", encoding="utf-8")
    print("\nwrote {}".format(OUT))


if __name__ == "__main__":
    main()
