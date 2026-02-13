# Kinase–Phosphosite Heterogeneous Graph

This project builds a heterogeneous network where nodes represent:
- kinases
- phosphatases
- proteins (substrates)
- phosphosites (protein-site)

Edges represent functional associations derived from multiple data sources:
- kinase/phosphatase targeting of sites
- protein-to-site containment
- site-site pathway co-membership (PTMsigDB)
- protein-protein pathway co-membership (MsigDB)
- protein-protein physical/functional interactions (high_confidence_score)
- site-site coevolution links (PTMcode2)

The output is a pair of files:
- `outputs/nodes.csv.gz`
- `outputs/edges.csv.gz`

These are designed to be loaded into Python, Neo4j, Cytoscape, Gephi, etc.

---

## Project Layout

kinase-phosphosite-graph/
data/ # input datasets (not recommended to push to GitHub)
outputs/ # generated outputs (optional to push)
src/
main.py # pipeline entrypoint
config.py # paths and constants
io_utils.py # file readers
normalize.py # normalization rules for proteins and sites
sanity_checks.py # basic validation checks
post_build_summary.py # post-build report on outputs
builders/
kinase_substrate.py
ppase_substrate.py
ptmsigdb.py
msigdb.py
ppi.py
ptmcode2.py

---

## Inputs

Place these files inside `data/`:

1. `Kinase_Substrate_Dataset`  
2. `PPase_protSubtrates_201903.xls`  
3. `PTMsigDB.txt`  
4. `MsigDB.txt`  
5. `high_confidence_score.csv` (or your equivalent)  
6. `gene_protein.csv` (Ensembl protein pairs -> gene symbols)  
7. `PTMcode2_associations_between_proteins.txt` (PTMcode2 v2, 2014 format)

Paths are controlled in `src/config.py`.

---

## Node Schema

`outputs/nodes.csv.gz` columns:
- `node_id` (string)
- `node_type` (string)

Node types used:
- `kinase`
- `phosphatase`
- `protein`
- `site`

### Site node ID rule
Phosphosite nodes are globally unique by construction:

`{PROTEIN}-{SITE_LABEL}`

Examples:
- `MAP2K1-S298`
- `FXYD1-S88`

This guarantees that `S298` in different proteins is not treated as the same node.

---

## Edge Schema

`outputs/edges.csv.gz` columns:
- `source` (string)
- `target` (string)
- `relation` (string)
- optional extra columns per builder (example: `confidence_score`)

Relations used:
- `has_site`
- `phosphorylates`
- `dephosphorylates`
- `site_same_pathway`
- `protein_same_pathway`
- `ppi_high_confidence`
- `site_coevolution`

### Meaning of key relations

- `PROTEIN -> SITE : has_site`  
  The site belongs to the protein.

- `KINASE -> SITE : phosphorylates`  
  Kinase is reported to phosphorylate that phosphosite.

- `PHOSPHATASE -> SITE : dephosphorylates`  
  Phosphatase is reported to dephosphorylate that phosphosite.

- `SITE -> SITE : site_same_pathway`  
  Two sites appear in the same PTMsigDB pathway set.

- `PROTEIN -> PROTEIN : protein_same_pathway`  
  Two proteins appear in the same MsigDB pathway set.

- `PROTEIN -> PROTEIN : ppi_high_confidence`  
  Protein-protein interaction edge derived from the high confidence PPI file after Ensembl->gene mapping.

- `SITE -> SITE : site_coevolution`  
  Site-site edge from PTMcode2 where:
  - Species is `Homo sapiens`
  - Both PTMs are `phosphorylation`
  - `Coevolution_evidence == 1`

---

## Build Rules Per Dataset

### 1) Kinase_Substrate_Dataset
Columns used (as in file):
- `KIN_ORGANISM`, `SUB_ORGANISM`
- `KINASE` (kinase name)
- `SUB_GENE` or substrate name (depending on file column naming)
- `SUB_MOD_RSD` (site)

Filters:
- keep rows where both kinase organism and substrate organism are human
- drop rows with missing kinase/substrate/site fields
- drop rows where multiple sites appear in one field (kept only rows with exactly one site token)

Nodes created:
- kinase node
- protein node (substrate)
- site node `PROTEIN-SITE`

Edges created:
- `PROTEIN -> SITE : has_site`
- `KINASE -> SITE : phosphorylates`

### 2) PPase_protSubtrates_201903
Columns used:
- phosphatase name
- substrate name
- site

Filters:
- drop rows with missing fields
- drop rows where multiple sites are reported
- normalize Ser/Thr/Tyr to S/T/Y

Nodes created:
- phosphatase node
- protein node
- site node `PROTEIN-SITE`

Edges created:
- `PROTEIN -> SITE : has_site`
- `PHOSPHATASE -> SITE : dephosphorylates`

### 3) PTMsigDB (site-level pathways)
Each line is: pathway name + a list of sites.

Rule:
- only keep sites that already exist from the base graph (Kinase + PPase)
- add clique edges among sites within the same pathway

Edge:
- `SITE -> SITE : site_same_pathway`

Safeguard:
- skip very large pathways to avoid edge explosion using `MAX_SITE_CLIQUE` (config)

### 4) MsigDB (protein-level pathways)
Each line is: pathway name + a comma-separated list of proteins.

Rule:
- only keep proteins already present in the base graph
- add clique edges among proteins within the same pathway

Edge:
- `PROTEIN -> PROTEIN : protein_same_pathway`

Safeguard:
- skip very large pathways using `MAX_PROTEIN_CLIQUE` (config)

### 5) High confidence PPI (Ensembl protein IDs)
File contains edges between Ensembl protein IDs.

Rule:
- map Ensembl protein IDs to gene symbols using `gene_protein.csv`
- if an Ensembl protein maps to multiple genes, we choose a single gene deterministically:
  - we currently pick the `max()` gene symbol (documented choice, can be revised later)
- only keep PPI edges where both proteins exist in our protein nodes

Edge:
- `PROTEIN -> PROTEIN : ppi_high_confidence`
Optional attribute:
- `confidence_score`

### 6) PTMcode2 associations
PTMcode2 file contains PTM pairs with evidence columns.

Filters:
- Species must be `Homo sapiens`
- PTM1 and PTM2 must be `phosphorylation`
- coevolution evidence must be `1`
- site must parse into a standard label (S/T/Y + number)
- site must exist in our site nodes

Edge:
- `SITE -> SITE : site_coevolution`

Note on file parsing:
- PTMcode2 has header/comment lines. The pipeline skips the header correctly.

---

## How to Run

Create and activate a virtual environment, then install dependencies:
pip install -r requirements.txt

```bash
python -m src.main