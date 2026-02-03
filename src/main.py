from config import *
import re
import gzip
import pandas as pd
from tqdm import tqdm


# -------------------------
# Helpers
# -------------------------
AA3_TO_AA1 = {"SER": "S", "THR": "T", "TYR": "Y"}

def norm_gene(x: str) -> str | None:
    if pd.isna(x):
        return None
    s = str(x).strip()
    return s.upper() if s else None

def norm_site(x: str) -> str | None:
    """
    Accepts: S298, Ser-88, Tyr-256, Thr-12 (with optional spaces)
    Returns: S298, Y256, T12
    """
    if pd.isna(x):
        return None
    s = str(x).strip()

    m = re.fullmatch(r"([STY])\s*[-]?\s*(\d+)", s, flags=re.I)
    if m:
        return m.group(1).upper() + m.group(2)

    m = re.fullmatch(r"(Ser|Thr|Tyr)\s*[-]?\s*(\d+)", s, flags=re.I)
    if m:
        aa = AA3_TO_AA1[m.group(1).upper()]
        return aa + m.group(2)

    return None

def clique_or_star(nodes: list[str], max_clique: int = 100):
    """
    Returns list of (a,b) undirected edges.
    If too big, returns star edges from first node.
    """
    if len(nodes) < 2:
        return []
    if len(nodes) > max_clique:
        hub = nodes[0]
        return [(hub, x) for x in nodes[1:]]
    out = []
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            out.append((nodes[i], nodes[j]))
    return out

# -------------------------
# Build edges
# -------------------------
def build_edges(
    kinase_tsv: str,
    phosphatase_xls: str,
    ptmsig_txt: str,
    msig_txt: str,
    ppi_csv: str,
    protein_map_csv: str,
    ptmcode_gz: str,
):
    edges = []

    # 1) Kinase-substrate (human-human only)
    ks = pd.read_csv(kinase_tsv, sep="\t")
    ks = ks[
        (ks["KIN_ORGANISM"].str.lower() == "human") &
        (ks["SUB_ORGANISM"].str.lower() == "human")
    ].copy()

    ks["KIN"] = ks["KINASE"].map(norm_gene)
    ks["SUB"] = ks["SUB_GENE"].map(norm_gene)
    ks["R"] = ks["SUB_MOD_RSD"].map(norm_site)
    ks = ks.dropna(subset=["KIN", "SUB", "R"])
    ks["SITE_ID"] = ks["SUB"] + "-" + ks["R"]

    for row in ks.itertuples(index=False):
        edges.append((row.SUB, row.SITE_ID, "has_site", "Kinase_Substrate_Dataset", None))
        edges.append((row.KIN, row.SITE_ID, "phosphorylates", "Kinase_Substrate_Dataset", None))

    # 2) Phosphatase-substrate
    # needs: pip install xlrd==2.0.1
    pp = pd.read_excel(phosphatase_xls, engine="xlrd")
    pp["PHOS"] = pp["Phosphatase entry names"].map(norm_gene)
    pp["SUB"] = pp["Substrate entry names"].map(norm_gene)
    pp["R"] = pp["Dephosphosites"].map(norm_site)

    # remove multi-site rows
    pp = pp[pp["R"].notna()]
    pp = pp[~pp["Dephosphosites"].astype(str).str.contains(r"[;,/]")]

    pp["SITE_ID"] = pp["SUB"] + "-" + pp["R"]
    pp = pp.dropna(subset=["PHOS", "SUB", "SITE_ID"])

    for row in pp.itertuples(index=False):
        edges.append((row.SUB, row.SITE_ID, "has_site", "PPase_protSubtrates", None))
        edges.append((row.PHOS, row.SITE_ID, "dephosphorylates", "PPase_protSubtrates", None))

    # Node universe is from kinase + phosphatase
    site_universe = set([e[1] for e in edges if e[2] in {"phosphorylates", "dephosphorylates"}])

    # 3) PTMsigDB site-level pathways
    with open(ptmsig_txt, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if ":" in line:
                pathway, rest = line.split(":", 1)
                pathway = pathway.strip()
                tokens = rest.strip().split()
            else:
                pathway = None
                tokens = line.split()

            sites = []
            for t in tokens:
                m = re.fullmatch(r"([A-Za-z0-9\-]+)-([STY])(\d+)", t)
                if not m:
                    continue
                gene = norm_gene(m.group(1))
                r = m.group(2).upper() + m.group(3)
                sid = f"{gene}-{r}"
                if sid in site_universe:
                    sites.append(sid)

            for a, b in clique_or_star(sites, max_clique=50):
                edges.append((a, b, "site_association_pathway", "PTMsigDB", pathway))

    # 4) MsigDB protein-level pathways
    protein_universe = set()
    for s, t, rel, _, _ in edges:
        if rel == "has_site":
            protein_universe.add(s)

    with open(msig_txt, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            line = line.strip()
            if not line or "\t" not in line:
                continue
            pathway, genes_str = line.split("\t", 1)
            genes = [norm_gene(g) for g in genes_str.split(",")]
            genes = [g for g in genes if g and g in protein_universe]

            for a, b in clique_or_star(genes, max_clique=100):
                edges.append((a, b, "protein_association_pathway", "MsigDB", pathway))

    # 5) PPI edges (Ensembl protein IDs -> gene symbols)
    # Build a protein_id -> gene symbol map from the large mapping file (chunked)
    ppi_ids = set()
    for chunk in pd.read_csv(ppi_csv, usecols=["protein1", "protein2"], chunksize=200000):
        ppi_ids.update(chunk["protein1"].unique().tolist())
        ppi_ids.update(chunk["protein2"].unique().tolist())

    prot_to_gene = {}
    for chunk in tqdm(pd.read_csv(protein_map_csv, chunksize=300000), desc="Mapping proteins"):
        m1 = chunk["protein1"].isin(ppi_ids)
        if m1.any():
            for p, g in zip(chunk.loc[m1, "protein1"], chunk.loc[m1, "gene1"]):
                prot_to_gene.setdefault(p, norm_gene(g))
        m2 = chunk["protein2"].isin(ppi_ids)
        if m2.any():
            for p, g in zip(chunk.loc[m2, "protein2"], chunk.loc[m2, "gene2"]):
                prot_to_gene.setdefault(p, norm_gene(g))

    for chunk in tqdm(pd.read_csv(ppi_csv, chunksize=300000), desc="PPI edges"):
        g1 = chunk["protein1"].map(prot_to_gene)
        g2 = chunk["protein2"].map(prot_to_gene)
        mask = g1.notna() & g2.notna()
        sub = chunk.loc[mask, ["confidence_score"]].copy()
        sub["g1"] = g1[mask]
        sub["g2"] = g2[mask]
        sub = sub[(sub["g1"].isin(protein_universe)) & (sub["g2"].isin(protein_universe))]

        for row in sub.itertuples(index=False):
            edges.append((row.g1, row.g2, "protein_interaction", "high_confidence_score", int(row.confidence_score)))

    # 6) PTMcode2 coevolution (optional, often limited overlap)
    with gzip.open(ptmcode_gz, "rt", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                continue
            prot1, prot2, species, ptm1, res1, *_rest = parts
            ptm2 = parts[7]
            res2 = parts[8]
            coevo = parts[11]

            if species != "Homo sapiens":
                continue
            if coevo != "1":
                continue
            if ptm1 != "phosphorylation" or ptm2 != "phosphorylation":
                continue

            g1 = norm_gene(prot1)
            g2 = norm_gene(prot2)
            r1 = norm_site(res1)
            r2 = norm_site(res2)
            if not (g1 and g2 and r1 and r2):
                continue

            s1 = f"{g1}-{r1}"
            s2 = f"{g2}-{r2}"
            if s1 in site_universe and s2 in site_universe and s1 != s2:
                edges.append((s1, s2, "site_association_coevolution", "PTMcode2", None))

    return pd.DataFrame(edges, columns=["source", "target", "relation", "evidence", "meta"])

# -------------------------
# Run + export
# -------------------------
if __name__ == "__main__":
    df_edges = build_edges(
        kinase_tsv="data/Kinase_Substrate_Dataset",
        phosphatase_xls="data/PPase_protSubtrates_201903.xls",
        ptmsig_txt="data/PTMsigDB.txt",
        msig_txt="data/MsigDB.txt",
        ppi_csv="data/high_confidence_score(in).csv",
        protein_map_csv="data/gene_protein(in).csv",
        ptmcode_gz="data/PTMcode2_associations_between_proteins.txt.gz",
    )

    df_edges.to_csv("outputs/edges.csv.gz", index=False, compression="gzip")
    print(df_edges["relation"].value_counts())
    print("Saved outputs/edges.csv.gz")
