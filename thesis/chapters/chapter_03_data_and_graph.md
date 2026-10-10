# CHAPTER III: DATA SOURCES AND GRAPH CONSTRUCTION

> Draft. Every count in this chapter was recomputed directly from the artefacts
> in `outputs/` rather than copied from project documentation; where the two
> disagree, the recomputed value is used and the discrepancy is noted.

## 3.1 Overview

The modelling substrate for this thesis is a heterogeneous graph whose nodes are
proteins and phosphorylation sites and whose edges encode several distinct kinds
of biological evidence. Construction proceeds in two stages. A **generic graph**
is assembled from curated public resources and represents context-free prior
knowledge. A **liver-extended graph** then adds proteins and sites observed in
liver cancer tissue together with correlation edges derived from liver
phosphoproteomic measurements. Four experimental networks are obtained from
these two files by selecting which edge types participate, so no graph is ever
rebuilt per experiment.

This chapter describes the inputs, the node and edge schema, the construction
procedure, the resulting statistics, and the properties of the constructed graph
that bear on the results in Chapter V.

## 3.2 Data Sources

Seven curated datasets supply the generic graph, and one proteomics dataset
supplies the liver context.

| Source | Contribution | File |
|---|---|---|
| PhosphoSitePlus | kinase-substrate edges; site definitions | `Kinase_Substrate_Dataset.tsv` |
| Phosphatase-substrate table | dephosphorylation edges | `PPase_protSubtrates_201903.xls` |
| PTMsigDB | site-level pathway co-membership | `PTMsigDB.txt` |
| MSigDB | protein-level pathway co-membership | `MsigDB.txt` |
| High-confidence PPI | protein-protein interactions | `high_confidence_score.csv` |
| PTMcode2 | site-site co-evolution | `PTMcode2_associations_between_proteins.txt.gz` |
| Gene-protein map | identifier normalization | `gene_protein.csv` |
| Liver cancer proteomics | liver nodes, correlation edges | `LiverCancer_ProtExp_Phospho_casecntrl.xlsx` |

The liver dataset provides matched protein expression and phosphorylation
measurements for **18 tumor** and **18 control** samples. It contributes two
sheets, `ProteinExpression` and `Phosphorylation`; only the phosphorylation
sheet is used to build correlation edges, for reasons given in Section 3.6.4.

Identifier normalization is applied before merging so that gene symbols and site
labels from different sources resolve to the same node. Large inputs are cached
in preprocessed form under `data/processed/` so the build is reproducible
without reparsing the raw files.

## 3.3 Node Schema

The graph uses exactly **two node types**:

- `PROTEIN:<GENE>` — any protein, including kinases and phosphatases
- `SITE:<GENE>-<residue><position>` — an individual phosphorylation site

Examples from the constructed graph are `PROTEIN:EIF2AK1`, `SITE:EIF2S1-S52`,
and `SITE:A1CF-T14`.

Enzymatic role is deliberately **not** encoded as a node type. Kinases and
phosphatases are proteins, and their roles are carried by node attributes
(`is_kinase`, `is_phosphatase`, `protein_role`) and by the type of the edges
they participate in. This choice keeps the topology simple while preserving
interpretability: every interaction involving a protein, whether catalytic or
not, contributes to a single unified representation of that protein rather than
being split across role-specific duplicates.

Modelling at site rather than protein granularity is essential, since the
prediction target is a kinase-to-site association and different sites on one
protein may be controlled by different kinases.

## 3.4 Edge Schema

Seven base relation types are shared by both graphs, and six correlation
relations exist only in the liver graph.

**Base relations**

| Relation | Endpoints | Meaning |
|---|---|---|
| `phosphorylates` | protein → site | kinase acts on site (**prediction target**) |
| `dephosphorylates` | protein → site | phosphatase acts on site |
| `has_site` | protein → site | site belongs to this protein |
| `protein_same_pathway` | protein ↔ protein | shared pathway membership |
| `ppi_high_confidence` | protein ↔ protein | high-confidence interaction |
| `site_same_pathway` | site ↔ site | shared pathway membership |
| `site_coevolution` | site ↔ site | co-evolution evidence |

**Liver correlation relations**, each split into positive and negative:
`site_corr_ctrl_pos/neg`, `site_corr_cancer_pos/neg`, `site_corr_fc_pos/neg`.

The `phosphorylates` relation requires emphasis because it plays two roles. It
is the supervised label, and it is also an edge in the graph. Chapter IV
describes how it is removed from the graph before any embedding is computed, to
prevent the answer leaking into the node representations.

Pathway co-membership is expressed as a clique over all members of a pathway,
which can produce very large edge sets. Clique size is therefore capped at 50
members for sites and 100 for proteins, preventing a single large pathway from
dominating the edge count.

## 3.5 Generic Graph Construction

The generic build loads PhosphoSitePlus and the phosphatase table, creates
protein and site nodes with role metadata, then layers on pathway
co-membership, protein interactions, and co-evolution edges. PhosphoSitePlus
rows are filtered to human entries with exactly one site token.

The resulting graph contains **13,318 nodes** (10,196 sites and 3,122 proteins)
and **1,095,542 edges**, with density 0.0124 and a single connected component.
Edges decompose as follows.

| Relation | Count |
|---|---|
| `protein_same_pathway` | 992,485 |
| `ppi_high_confidence` | 36,009 |
| `site_coevolution` | 24,005 |
| `site_same_pathway` | 17,922 |
| `phosphorylates` | 14,596 |
| `has_site` | 10,196 |
| `dephosphorylates` | 329 |

Pathway co-membership supplies 90.6% of all edges, a concentration that matters
for interpreting any result driven by graph connectivity.

> **Discrepancy noted.** Earlier project documentation and the conference poster
> report 1,070,352 generic edges and density 0.0121. The current
> `outputs/edges.csv.gz` contains 1,095,542 edges summing exactly to the table
> above. The recomputed figures are used throughout this thesis.

## 3.6 Liver Context Integration

### 3.6.1 Adding Liver-Observed Entities

Proteins and phosphosites detected in the liver dataset but absent from the
generic graph are added as new nodes, and the pathway, interaction, and
co-evolution integration steps are re-run over the expanded node set so that new
entities are connected by the same evidence types as existing ones. This grows
the graph from 13,318 to **16,605 nodes** (13,410 sites and 3,195 proteins). The
increase is overwhelmingly on the site side, reflecting the depth of
phosphoproteomic measurement relative to the curated annotation.

The base (non-correlation) edge count rises correspondingly to **1,133,920**.

### 3.6.2 Correlation Edge Construction

Correlation edges connect pairs of phosphosites whose measured abundances move
together across samples. For each pair of sites, a Pearson correlation is
computed over the samples where both are observed, requiring a minimum of six
overlapping measurements. Correlations above the **80th percentile** of the
positive distribution become positive edges, and correlations below the
corresponding negative threshold become negative edges. Positive and negative
relations are kept distinct rather than merged, so that anti-correlated and
co-correlated structure remain separable.

An 85th-percentile threshold was tested on the fold-change variant during
development and performed worse, so the 80th percentile is retained for all
three correlation types.

### 3.6.3 Three Correlation Types

Three correlation edge sets are computed in a single pass and stored together in
one file:

- **Control** — raw abundance correlation across the 18 healthy samples,
  capturing baseline co-regulation in non-tumor liver tissue.
- **Cancer** — raw abundance correlation across the 18 tumor samples,
  capturing co-regulation in the disease state.
- **Liver FC** — correlation of tumor-versus-control fold-change profiles,
  capturing how site regulation *changes* in tumor relative to baseline.

Storing all three in one artefact and selecting among them at run time is a
deliberate design decision: it guarantees that the four experimental networks
differ only in which correlation edges are active, with the base graph held
exactly constant.

### 3.6.4 Protein Co-Expression Edges Are Disabled

The construction code contains logic for protein-level correlation edges, but it
is switched off. After preprocessing, only three usable protein rows remain in
the input, which is far too few to estimate a meaningful correlation structure.
No protein correlation edges appear in any network used in this thesis.

### 3.6.5 Sample-Size Characteristics of the Correlation Edges

Because the correlation is computed over the samples where both sites are
observed, and phosphoproteomic data are sparse, the effective sample size varies
widely across pairs. For the 976,855 Liver FC correlation edges:

| Overlapping samples | Edges | Share |
|---|---|---|
| 6–7 | 256,805 | 26.3% |
| 8–11 | 340,267 | 34.8% |
| 12–14 | 174,904 | 17.9% |
| 15–17 | 126,741 | 13.0% |
| 18 (complete) | 78,138 | **8.0%** |

Only 8.0% of correlation edges rest on the full complement of 18 samples, and
the most common value is the minimum permitted, n = 6 (151,197 edges). A Pearson
correlation estimated from six observations has very high variance, so a
substantial share of these edges may reflect noise rather than co-regulation.

It is worth stating what this does **not** show. Mean |r| among retained edges
is essentially flat across sample-size buckets and is in fact highest at n = 18
(0.856 at n = 6–7 against 0.875 at n = 18), so the percentile threshold does not
preferentially admit low-sample pairs. The concern is estimator variance, not
systematic inflation. This measurement is conditional on edges that already
passed the 80th-percentile cut; the unconditional distribution was not examined.

## 3.7 The Four-Network Design

Four networks are defined by edge-type selection over the two graph files.

| Network | File | Correlation edges active |
|---|---|---|
| Generic | generic | none |
| Control | liver | `site_corr_ctrl_pos/neg` |
| Cancer | liver | `site_corr_cancer_pos/neg` |
| Liver FC | liver | `site_corr_fc_pos/neg` |

Because the embedding is recomputed independently for each network from that
network's own graph, context is encoded in the geometry of the representation
rather than appended as a post-hoc feature.

## 3.8 Graph Statistics

| Property | Generic | Control | Cancer | Liver FC |
|---|---|---|---|---|
| Nodes | 13,318 | 16,605 | 16,605 | 16,605 |
| Sites / proteins | 10,196 / 3,122 | 13,410 / 3,195 | 13,410 / 3,195 | 13,410 / 3,195 |
| Base edges | 1,095,542 | 1,133,920 | 1,133,920 | 1,133,920 |
| Correlation edges | 0 | 961,364 | 1,015,615 | 976,855 |
| **Total edges** | **1,095,542** | **2,095,284** | **2,149,535** | **2,110,775** |
| Density | 0.0124 | 0.0152 | 0.0156 | 0.0153 |
| Connected components | 1 | 1 | 1 | 1 |

Each network forms a single connected component, confirming that the
multi-layer integration succeeds in linking the curated and measured evidence
rather than producing disjoint islands.

> **Reporting note.** Project documentation describes the liver graph as having
> "2.1 million edges". That figure corresponds to the Liver FC variant
> specifically (2,110,775). The liver edge *file* contains 4,087,754 rows,
> because it stores all three correlation types simultaneously. Per-variant
> counts are used throughout this thesis.

## 3.9 Annotation Coverage of the Correlation Edges

A property of the correlation layer that becomes important in Chapter V is how
its edges relate to annotated biology. Of 10,146 sites carrying at least one
PhosphoSitePlus kinase, 9,258 are substrates of a well-annotated kinase, 1,230
of a moderately annotated kinase, and only 360 of a poorly annotated kinase.

Classifying the Liver FC correlation edges by their endpoints:

| Endpoint property | Edges | Share |
|---|---|---|
| Touches a well-annotated kinase's substrate | 193,271 | 19.8% |
| Touches a poorly annotated kinase's substrate | 12,079 | 1.2% |
| Neither endpoint has any kinase annotation | 768,829 | **78.7%** |

Nearly four-fifths of the correlation edges connect sites with no kinase
annotation at all, and among annotated endpoints the ratio of well-annotated to
poorly annotated exposure is roughly sixteen to one. The context layer therefore
adds most of its connectivity in a region of the graph where no supervision
exists, and what annotated exposure it does add is concentrated on kinases that
were already data-rich.

## 3.10 Summary

The construction yields a single connected heterogeneous graph of 16,605 nodes
in which curated kinase-substrate, phosphatase, pathway, interaction, and
co-evolution evidence is combined with liver cancer phosphoproteomic
correlations, and from which four controlled network variants can be selected
without rebuilding. Two properties established here constrain the
interpretation of later results: the correlation layer is built largely from
small-sample estimates, and it attaches predominantly to unannotated regions of
the graph. Chapter IV describes how predictions are made and evaluated on these
networks.
