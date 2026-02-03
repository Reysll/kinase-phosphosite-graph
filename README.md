# kinase-phosphosite-graph

## Overview
This project constructs a heterogeneous graph for kinase–substrate–phosphosite relationships in human phosphorylation data.

The goal is to generate a clean node/edge representation that can be used for downstream analysis and machine learning (e.g., graph-based prediction tasks).

## Current status
Implemented:
- Graph construction from `Kinase_Substrate_Dataset` (human-only filtering)
- Node and edge export to compressed CSV
- Sanity checks to ensure phosphosites are uniquely identified by protein context

Planned (later):
- Add pathway edges from PTMsigDB (site-level pathways)
- Add pathway edges from MSigDB (protein-level pathways)
- Add PPI edges using `high_confidence_score(in).csv` with Ensembl-to-gene mapping
- Add site coevolution edges from PTMcode2 (Ensembl-to-gene mapping required)

## Graph schema

### Node types
Nodes are typed and stored in `outputs/nodes.csv.gz` with columns:
- `node_id`
- `node_type`

| node_type | meaning |
|----------|---------|
| kinase   | kinase enzyme |
| protein  | substrate protein |
| site     | phosphosite located on a protein |

### Phosphosite node ID rule (important)
Phosphosite labels like `S298` are NOT globally unique.

Therefore, site nodes are uniquely named as:
`<PROTEIN>-<AA><POSITION>`

Examples:
- `MAP2K1-S298`
- `FXYD1-S88`

This prevents collisions where two different proteins share the same site label.

### Edge types
Edges are stored in `outputs/edges.csv.gz` with columns:
- `source`
- `target`
- `relation`

| relation | source -> target | meaning |
|----------|------------------|---------|
| has_site | protein -> site  | site belongs to a protein |
| phosphorylates | kinase -> site | kinase phosphorylates that site |

## Data sources
Raw data files are stored locally in `data/` and are NOT committed to GitHub.

Current input used:
- `Kinase_Substrate_Dataset` (filtered to human kinase and human substrate rows)

## Outputs
Generated files:
- `outputs/nodes.csv.gz`
- `outputs/edges.csv.gz`

## Sanity checks
We run checks after building the graph to verify:
1. Site nodes follow the `PROTEIN-S###` naming rule
2. The same site label (example: `S298`) can appear on multiple proteins (expected)
3. No duplicate sites exist within the same protein context

## How to run

Create/activate your virtual environment, install requirements, then:

```bash
python -m src.main
