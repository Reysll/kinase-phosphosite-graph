# Graph construction and integration: a step-by-step walkthrough

Prepared for a meeting with Dr. Ayati. Written to be read without opening the
code, but every claim names the file and function that produces it so anything
can be checked on the spot.

**Everything in this document was recomputed read-only on 2026-10-25** from the
raw files in `data/` and the built graph in `outputs/`. Nothing was rebuilt, and
no code, data or result file was modified. Where a number could only be
recovered by rerunning the build (which would overwrite `outputs/`), that is
stated instead of estimated.

Interpreter for every command: `.venv312/Scripts/python`.

---

## 1. What the graph is

The graph has exactly **two kinds of node**.

| Node type | Written as | Example | Meaning |
|---|---|---|---|
| Protein | `PROTEIN:<GENE>` | `PROTEIN:AKT1` | a protein, named by its gene symbol |
| Phosphosite | `SITE:<GENE>-<residue><position>` | `SITE:EIF2S1-S52` | one specific phosphorylatable residue on one protein |

A site belongs to exactly one protein. `SITE:EIF2S1-S52` is serine 52 of EIF2S1.

There are **seven kinds of edge**, each coming from one source dataset.

| Edge type | Connects | Biological meaning | Source |
|---|---|---|---|
| `has_site` | protein to site | "this residue is on this protein" | PhosphoSitePlus, DEPOD |
| `phosphorylates` | protein to site | "this kinase is known to phosphorylate this site" | PhosphoSitePlus |
| `dephosphorylates` | protein to site | "this phosphatase removes this phosphate" | DEPOD |
| `ppi_high_confidence` | protein to protein | "these two proteins physically interact" | PPI file |
| `protein_same_pathway` | protein to protein | "both proteins are in the same annotated gene set" | MSigDB |
| `site_same_pathway` | site to site | "both sites are in the same signature set" | PTMsigDB |
| `site_coevolution` | site to site | "these two sites co-evolve across species" | PTMcode2 |

**The `phosphorylates` edges are the prediction target.** They are deliberately
removed before the embedding is computed, so the model never sees a
kinase-substrate link while learning the structure, and they are used only as
training labels. This is the leakage guard, in
`src_prediction/engine/leave_one_out.py`.

---

## 2. The single most important structural fact

Only **two** of the seven sources create nodes. The other five can only add
edges between nodes that already exist.

Look at the build order in `src/runners/main.py:100-190`:

1. PhosphoSitePlus is parsed. It creates protein nodes and site nodes.
2. DEPOD is parsed. It creates more protein nodes and site nodes.
3. The two are merged. **The node set is now final.**
4. PTMsigDB is called as `build_ptmsigdb_edges(PTMSIGDB_FILE, site_ids)`.
5. MSigDB is called as `build_msigdb_edges(MSIGDB_FILE, protein_ids)`.
6. PPI is called with `existing_protein_ids=protein_ids`.
7. PTMcode2 is called with `existing_site_ids=site_ids`.

Steps 4 to 7 each receive the already-built node set and discard any record
whose endpoints are not in it. Their return values are merged into the edge
table only, never into the node table.

**Consequence.** The universe of proteins and sites the model can ever reason
about is defined entirely by PhosphoSitePlus plus DEPOD. If a liver-relevant
site is absent from PhosphoSitePlus, no amount of pathway or interaction data
will bring it into the generic graph. This matters for the collaborator's
validation request, because the substrates they want are measured in the liver
proteomics, not necessarily annotated in PhosphoSitePlus.

---

## 3. Source by source

### 3.1 PhosphoSitePlus, `data/raw/Kinase_Substrate_Dataset.tsv`

The kinase-substrate annotation backbone.

| Stage | Rows | Where | Why |
|---|---|---|---|
| Raw file | **23,820** | | all species |
| Human kinase and human substrate | **14,676** | `src/graph/preprocess.py`, `preprocess_kinase_substrate()` | we model human signalling only; cached to `data/processed/Kinase_Substrate_Dataset_human.tsv` |
| Species check repeated | no change | `src/builders/kinase_substrate.py` | redundant after the step above, harmless |
| Exactly one residue token in the site field | not recoverable | `_has_exactly_one_site_token()` | a row naming two residues cannot be turned into one site node, so it is dropped rather than guessed |
| Name and site must normalize | not recoverable | `normalize_protein()`, `normalize_site()` | discards blank or unparseable identifiers |

**Contributed to the final graph: 14,596 `phosphorylates` edges and 10,196
`has_site` edges**, from 420 distinct kinases.

The two "not recoverable" rows are honest gaps. The builder prints its own drop
counts when it runs, but those logs were not saved, and recovering them would
require rerunning the build, which overwrites `outputs/`. The endpoints are
solid: 14,676 human rows in, 14,596 unique kinase-to-site edges out.

### 3.2 DEPOD, `data/raw/PPase_protSubtrates_201903.xls`

This is the source of the discrepancy Dr. Ayati found, so it is traced in full.

| Stage | Rows | Distinct phosphatases | Where | Why |
|---|---|---|---|---|
| Raw file | **938** | **115** | | DEPOD protein-substrate table, March 2019 release |
| Phosphatase, substrate and dephosphosite all present | **517** | **82** | `src/builders/ppase_substrate.py` | many DEPOD records name a substrate protein but no specific residue, so no site node can be made |
| Exactly one residue token, then deduplicated | **329 edges** | **65** | `_has_exactly_one_site_token()` | rows listing several dephosphosites are discarded whole |

**Contributed to the final graph: 329 `dephosphorylates` edges from 65
phosphatases onto 171 sites.**

Three things follow.

1. **Dr. Ayati's count is correct.** The file contains 115 distinct
   phosphatases, which matches her 114 to within one entry. The 250 reported
   earlier is not reproducible from this file. The most likely explanation is
   that 250 counts DEPOD's full phosphatase catalogue, most of which have no
   annotated protein substrate and so never appear in this table, but that is a
   hypothesis and should be checked against whichever DEPOD file was used then.
2. **Our pipeline is lossier than the file.** We end with 65, not 115. The
   larger loss is step two, DEPOD records without a named residue. The
   single-site filter then removes a further 17 phosphatases.
3. **The blast radius is very small.** 329 edges out of 1,095,542 is **0.03%**
   of the generic graph. Fixing the phosphatase layer is correct housekeeping,
   but it cannot explain or change any result we have measured. This is worth
   saying plainly, because the DEPOD discrepancy is a good reason to audit the
   preprocessing and a poor reason to distrust the modelling results.

An independent check confirms the number. The graph stores a `protein_role`
column on every protein node, written by `_add_protein_metadata()` in
`src/runners/main.py:62`. Counting it in `outputs/nodes.csv.gz` gives 64
phosphatase plus 1 protein that is both a kinase and a phosphatase, so 65, and
419 kinase plus that same 1, so 420. Both match the edge counts exactly.

### 3.3 MSigDB, `data/raw/MsigDB.txt`

**This is 90.6% of the graph and deserves the most scrutiny.**

| Stage | Count | Where | Why |
|---|---|---|---|
| Raw file | **4,761** gene sets | | one pathway or signature per line |
| Restricted to proteins already in the graph | | `build_msigdb_edges(path, allowed_protein_ids)` | cannot connect a protein the graph does not have |
| Gene sets with more than 100 surviving members skipped | | `MAX_PROTEIN_CLIQUE = 100` in `src/core/config.py` | a 500-member set would alone add 124,750 edges |

**Contributed: 992,485 `protein_same_pathway` edges.**

The mechanism is the thing to understand. For each gene set, every surviving
member is connected to every other member, a **clique**. A set with 100 members
contributes 4,950 edges by itself. This is why one source accounts for nine
tenths of the graph.

Here is the whole generic graph by edge type:

| Relation | Edges | Share |
|---|---|---|
| `protein_same_pathway` | 992,485 | **90.6%** |
| `ppi_high_confidence` | 36,009 | 3.3% |
| `site_coevolution` | 24,005 | 2.2% |
| `site_same_pathway` | 17,922 | 1.6% |
| `phosphorylates` (the target) | 14,596 | 1.3% |
| `has_site` | 10,196 | 0.9% |
| `dephosphorylates` | 329 | 0.03% |
| **Total** | **1,095,542** | |

**Why this matters more than DEPOD.** Pathway co-membership is not a
measurement. It is an annotation, it is transitive by construction, and it is
densest for well-studied proteins, which are exactly the kinases that already
have the most substrates. A graph that is 90.6% pathway cliques is a graph in
which a node's connectivity largely reflects how well studied it is. That is
consistent with our finding that a degree-preserving randomization reproduces
most of our results, and with kinase identity alone explaining 99.97% of the
variance in predicted rank.

### 3.4 PTMsigDB, `data/raw/PTMsigDB.txt`

| Stage | Count | Where | Why |
|---|---|---|---|
| Raw file | **687** signature sets | | site-level signatures |
| Restricted to sites already in the graph | | `build_ptmsigdb_edges(path, allowed_site_ids)` | |
| Sets with more than 50 surviving sites skipped | | `MAX_SITE_CLIQUE = 50` | same clique explosion guard |

**Contributed: 17,922 `site_same_pathway` edges**, also as cliques.

### 3.5 PPI, `data/raw/high_confidence_score.csv` plus the identifier map

| Stage | Count | Where | Why |
|---|---|---|---|
| Raw file | **472,000** protein pairs | | |
| Confidence threshold | **no filtering applied** | `min_confidence_score=0` at `src/runners/main.py` | see the warning below |
| Ensembl protein ID mapped to gene symbol | | `_build_ensp_to_gene_map()` in `src/builders/ppi.py` | the PPI file uses `9606.ENSP...` IDs, the graph uses gene symbols |
| Both ends must already be in the graph | | `existing_protein_ids` | |

**Contributed: 36,009 `ppi_high_confidence` edges.**

Two findings here, both new.

**Finding 1: the "high confidence" filter is not ours.** The relation is named
`ppi_high_confidence`, but the code passes `min_confidence_score=0`, so our
pipeline applies no confidence cut at all. Whatever confidence selection exists
is whatever was already baked into the file before it reached us. The name
therefore describes an assumption, not a step we perform.

**Finding 2: the identifier map appears to have been truncated by Excel.**
`data/raw/gene_protein.csv` has exactly **1,048,576 lines**, which is exactly
2^20, which is exactly Excel's maximum worksheet row count. A file that happens
to land on that number is almost always a file that was opened in Excel and
saved back.

The consequences are measurable:

- The map covers **14,521** distinct Ensembl protein IDs.
- The PPI file references **16,185** distinct Ensembl protein IDs.
- **1,664 of them (10.3%) have no gene mapping at all**, so every PPI edge
  touching them is silently discarded at `rows_dropped_no_map`.
- Those unmapped IDs appear in **6,960 of 472,000 PPI rows (1.5%)**.

So the practical cost is modest, 1.5% of candidate interactions, because the
unmapped proteins are mostly low-degree ones. But the cause is a data-handling
accident rather than a scientific choice, the loss is invisible in any output we
have, and the same map is also passed to PTMcode2. This is the clearest concrete
support for Dr. Ayati's instinct that the preprocessing needs an audit.

### 3.6 PTMcode2, `data/raw/PTMcode2_associations_between_proteins.txt.gz`

The largest raw input by far.

| Stage | Rows | Where | Why |
|---|---|---|---|
| Raw file | **18,727,764** | | all species, all PTM type pairs |
| Human, phosphorylation on both ends, coevolution evidence present | **3,861,786** | `preprocess_ptmcode2()` in `src/graph/preprocess.py` | cached to `data/processed/ptmcode2_filtered_66eef5e0ed.tsv.gz` |
| Residue labels must parse to a single site each | | `_extract_phospho_site_label()` | |
| Both sites must already be in the graph | | `existing_site_ids` | this is the dominant filter, since PTMcode2 covers far more sites than PhosphoSitePlus annotates |

**Contributed: 24,005 `site_coevolution` edges.**

Two inconsistencies between the preprocessing step and the builder, found while
tracing this:

1. **The first four rows of the cached file are skipped.**
   `preprocess_ptmcode2()` writes the cached file with `header=False`, so it has
   no header and no comment lines. `build_ptmcode2_edges()` then reads it with
   `skiprows=4`. That argument was written for the raw file, which does begin
   with comment lines. Applied to the preprocessed file it discards four genuine
   data rows. The numerical effect is at most four edges, so it changes nothing,
   but the preprocess-to-builder contract is mismatched and the next person to
   touch either side needs to know.
2. **The coevolution filter is stricter in the builder.** Preprocessing keeps
   `Coevolution_evidence >= 1`; the builder keeps only `== 1`. These agree only
   if the column is strictly 0 or 1. If it is ever greater than 1, those records
   survive preprocessing and are then dropped. Worth one check.

### 3.7 `data/raw/gene_protein.csv`

Not a biological source. It is an identifier bridge, 1,048,575 data rows of
`protein1, protein2, gene1, gene2`, used by PPI and PTMcode2 to turn Ensembl
protein IDs into gene symbols. See the truncation finding in 3.5.

---

## 4. The generic graph, assembled

| Quantity | Value |
|---|---|
| Nodes | **13,318** |
| of which phosphosites | 10,196 |
| of which proteins | 3,122 |
| Edges | **1,095,542** |
| Density | 0.0124 |
| Connected components | 1 |

Protein nodes by role, from `outputs/nodes.csv.gz`:

| Role | Count |
|---|---|
| protein (neither kinase nor phosphatase) | 2,638 |
| kinase | 419 |
| phosphatase | 64 |
| both kinase and phosphatase | 1 |

---

## 5. The liver extension

The liver graph starts from the generic graph and adds sites measured in the
liver proteomics plus a large set of **correlation edges**.

### 5.1 The data file

`data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx`, two sheets, read with
`header=1` in `src/graph/liver_network.py:576`.

| Sheet | Original table caption inside the file | Control columns | Tumour columns |
|---|---|---|---|
| `ProteinExpression` | "Table S4 Detailed information of unique proteins identified in HCC tumour tissues" | 18 | 18 |
| `Phosphorylation` | "Table S11 Identified and quantified phosphopeptides, phosphosites and corresponding phosphoproteins in HCC tumour and adjacent non-tumour tissues" | 18 | 18 |

Columns are identified purely by their names ending in `Control` or `Sample`,
in `_find_control_sample_columns()`.

**Provenance is still unresolved, and this is a question for the meeting.** The
two captions say these are supplementary Tables S4 and S11 of some HCC
proteogenomics paper, but the paper is not recorded anywhere in the repository.
Gao et al. (2019, *Cell*) is ruled out on sample count: that study reports 159
patients, and this file has 18 and 18. **We should not cite a source for this
dataset until someone can name the paper these two tables came from.**

### 5.2 Are the 18 and 18 matched pairs?

This matters for the collaborator's differential analysis, so it was checked
directly. The phosphorylation sheet's columns are, in file order:

```
Control:  F1  F2  F3   F7  F8  F9   F13 F14 F15   F19 F20 F21   F31 F32 F33   F43 F44 F45
Tumour:   F4  F5  F6   F10 F11 F12  F16 F17 F18   F22 F23 F24   F34 F35 F36   F46 F47 F48
```

The pattern is six blocks of three controls followed by three tumours. That is
the signature of a multiplexed run with six batches of six channels, three
tumour and three adjacent-normal per batch, which is consistent with **18
patients each contributing a tumour and its own adjacent non-tumour tissue**.
The sheet caption, "tumour and adjacent non-tumour tissues", says the same.

**But the pairing is not recorded.** Nothing in the file says which control
channel belongs to the same patient as which tumour channel. Within a batch,
`F1` could pair with `F4`, `F5` or `F6`. The protein sheet uses an entirely
different numbering (`F26` to `F61`), so the two sheets cannot be cross-walked
either.

So: the design is almost certainly paired, and a paired test is almost certainly
the right one, but **we cannot currently execute it**, because we do not know the
pairs. Resolving this needs either the source publication's sample table or the
original submitter. Section 7 of the differential-analysis plan proposes an
empirical fallback.

### 5.3 How correlation edges are built, and the thing Dr. Ayati asked about

Dr. Ayati asked: *"Don't we calculate correlation across all samples?"*

**The answer is no, and this is the most important technical point in this
document.**

Here is the whole mechanism, in plain terms.

A phosphoproteomics run does not detect every site in every sample. The matrix
is full of holes. Our pipeline handles that in two steps.

**Step one, which sites to keep.** A site is kept if it was quantified in at
least **6** of the 18 samples. That is `min_identified = 6`.

**Step two, how to correlate two kept sites.** The single line that does it is
`src/graph/liver_network.py:323`:

```python
corr = fc_mat.T.corr(method="pearson", min_periods=min_common)
```

pandas computes this using **pairwise complete observations**. For each *pair*
of sites it uses only the samples where **both** sites were measured, and
requires at least `min_common = 6` such overlapping samples. Neither value is
overridden anywhere: `src/runners/run_liver.py` passes only
`site_corr_percentile=80.0`, so both stay at their default of 6.

So two sites each measured in 8 samples, overlapping in 6, get a Pearson
correlation computed from **6 points**, not 18.

Here is how often that happens, counted directly from the 976,855 fold-change
correlation edges actually stored in `outputs/Liver_network_edges.csv.gz`:

| Overlapping samples | Edges | Share |
|---|---|---|
| 6 | 151,197 | **15.5%** (the most common value) |
| 7 | 105,608 | 10.8% |
| 8 | 89,868 | 9.2% |
| 9 | 103,346 | 10.6% |
| 10 | 76,937 | 7.9% |
| 11 | 70,116 | 7.2% |
| 12 | 71,173 | 7.3% |
| 13 | 52,233 | 5.3% |
| 14 | 51,498 | 5.3% |
| 15 | 61,914 | 6.3% |
| 16 | 27,808 | 2.8% |
| 17 | 37,019 | 3.8% |
| 18 (complete) | 78,138 | **8.0%** |

Grouped: 26.3% of edges rest on 6 or 7 samples, 61.1% on fewer than 12, and only
8.0% on all 18.

**Why this is a problem, stated carefully.** It is not that low-sample edges are
preferentially selected. We checked that and it is false: the mean absolute
correlation is essentially flat across sample-size buckets and is in fact
highest at n = 18. The problem is different and more basic:

The threshold for keeping an edge is the **80th percentile of all correlations
pooled together**, computed in `_build_corr_edges_from_fc_matrix()` with no
reference to how many samples each one used. A correlation of 0.9 from 6 points
and a correlation of 0.9 from 18 points are treated as equally trustworthy and
compete for the same cut. They are not equally trustworthy. With 6 points,
r = 0.9 is roughly p = 0.01; with 18 points it is roughly p = 1e-7. **No
significance test or multiple-testing correction is applied at any stage of edge
selection.**

The consequence is that an unknown but large share of the 976,855 correlation
edges are plausibly noise, and because the model's behaviour turns out to depend
mainly on how many edges a node gains rather than which edges, that noise is not
harmless.

### 5.4 The three correlation flavours and the four networks

One liver file holds all three correlation types at once, 4,087,754 edge rows.
Experiments select among them with `allowed_relations`.

| Network | Correlation edges used | Correlation edges | Total edges |
|---|---|---|---|
| Generic | none | 0 | 1,095,542 |
| Control | `site_corr_ctrl_pos/neg`, the 18 adjacent-normal samples | 961,364 | 2,095,284 |
| Cancer | `site_corr_cancer_pos/neg`, the 18 tumour samples | 1,015,615 | 2,149,535 |
| Liver FC | `site_corr_fc_pos/neg`, tumour divided by mean control | 976,855 | 2,110,775 |

Liver graph node counts: **16,605** nodes, 13,410 sites and 3,195 proteins.

Note that the liver base is 1,133,920 edges, not the generic 1,095,542, because
adding liver-measured sites also adds their `has_site` edges.

One more fact worth having in the room: **77% to 79% of the correlation edges
touch no PhosphoSitePlus-annotated site at either end.** They connect measured
liver sites to each other, in a region of the graph where no kinase-substrate
label exists. From `thesis/analysis/out/t13_endpoint_categories.csv`.

### 5.5 Protein co-expression edges are off

`liver_network.py` can also build protein-level correlation edges. It is
disabled, because the input had only three usable rows. Do not re-enable without
new data.

---

## 6. Discrepancies between the code and the documentation

Found while tracing. The first five were corrected in the docs on 2026-10-10
(commit `c5f6111`); the last three are new in this document.

| # | Claim in docs | What the code and data say | Status |
|---|---|---|---|
| 1 | Site IDs look like `SITE:GENE_pS473` | Actual format is `SITE:GENE-S473`, hyphen, no `p` | fixed |
| 2 | "Top-10% / Top-20%" metrics | These are Hit@10 and Hit@20, absolute rank cutoffs, not percentiles | fixed |
| 3 | 420 candidate kinases scored | 420 defined, **415** actually scored, 5 lack an embedding | fixed |
| 4 | KSEA gives 17 of 64 significant | 17 used a one-tailed p; two-tailed gives **12 of 64** | fixed |
| 5 | 6,909 LOO trials | 6,909 frozen, **15 unrankable in every run**, so all metrics use n = **6,894** | fixed |
| 6 | DEPOD contributes the phosphatase layer | 115 phosphatases in the file, **65** reach the graph | **new** |
| 7 | `ppi_high_confidence` implies a confidence filter | `min_confidence_score=0`, we apply none | **new** |
| 8 | identifier map is complete | exactly 2^20 lines, **1,664 of 16,185** PPI proteins unmapped | **new** |

---

## 6.1 A method defect outside the graph build: negative sampling

The audit that produced this document was extended to the prediction code, and
found one defect there. It is recorded here because it is the same class of
problem, a stated intent that the implementation does not carry out.

**What the code does.** Training the classifier requires negative examples,
meaning kinase and site pairs assumed not to interact. For each site the
selection is:

```python
negatives = [k for k in candidate_kinase_ids if k not in true_kinases_for_site]
if max_negatives is not None:
    negatives = negatives[:max_negatives]
```

with the cap set to 50. This takes the **first 50** entries of the candidate
list, not a random sample of 50. Whether that matters depends on the list's
order, and the candidate list is built with `.sort_values("node_id")`, so it is
**sorted alphabetically by gene symbol**: AAK1, AATK, ABL1, ABL2, ACVR1B, and so
on.

**Measured consequences**, reproducing the selection exactly over all 2,459
training sites:

| Quantity | Value |
|---|---|
| Distinct kinases ever used as a negative | **55 of 420** |
| Candidate-list positions they occupy | indices 0 to 54, the alphabetical head |
| Negative pairs generated in total | 122,950 |
| Times a typical such kinase is used | 2,458 or 2,459, i.e. nearly every site |
| Scored kinases that are **never** a negative | **361 of 415 (87.0%)** |

The classifier is therefore fitted against an alphabetically determined slice of
the candidate space. It never sees an explicit negative example for 87% of the
kinases it is later asked to rank.

**What has been tested, and what that test can and cannot see.** The
comparison below asks whether the 54 ever-negative kinases are pushed down
relative to the 361 never-negative ones. That is a **local** question. It
cannot answer the **global** one, which is more important: if the classifier
only ever learns to separate positives from the same 55 early-alphabet kinases,
it fits a direction in feature space that moves every kinase's score, including
those of kinases that never appear as negatives. A uniform shift of that kind
is invisible to a comparison between the two groups, because it moves both.
It is also a plausible mechanism for the near-fixed global kinase ordering
reported elsewhere in this work.

| Source | Ever-negative | Never-negative | p | Rank-biserial |
|---|---|---|---|---|
| Predict-all, Control | 213.5 | 207.0 | 0.740 | -0.028 |
| Predict-all, Liver FC | 215.5 | 207.0 | 0.698 | -0.033 |
| LOO Generic, adjusted rank | 240.0 | 191.5 | 0.060 | -0.170 |
| LOO Generic, within-category | 77.0 | 65.5 | 0.085 | -0.156 |
| LOO Control, adjusted rank | 206.0 | 196.8 | 0.569 | -0.052 |
| LOO Liver FC, adjusted rank | 206.0 | 197.3 | 0.563 | -0.052 |

No comparison reaches significance at 0.05. The Generic network shows the
largest gap, 48.5 rank positions in favour of never-negative kinases, but at
p = 0.060 across eight uncorrected comparisons it is not evidence of an effect.

**How to report it.** The defect is real and is a limitation of the training
design. On its consequences, claim **neither** bias **nor** absence of bias:
the local effect was not detected, and the global effect is untested. The only
way to test the global effect is to rerun with a uniformly random sampler and
compare. Such a rerun holds the frozen trial set, the candidate set, the
embeddings and the classifier fixed and changes only the negative selection, so
it is a valid paired comparison rather than a break in comparability. An
earlier version of this section treated the undetected local effect as evidence
of no effect, which overstated what the test could see.

Source: `thesis/analysis/out/negative_sampling_audit.txt`, reproduced by
`python -m thesis.analysis.negative_sampling_audit`.

---

## 7. What I would ask in the meeting

Ordered by how much they change what we do next.

1. **The dataset's provenance.** Which paper do Tables S4 and S11 come from? We
   cannot cite the liver data until this is settled, and Gao et al. is ruled out
   on sample count.
2. **The patient pairing.** Can we get the sample-to-patient table? The
   collaborator's differential analysis is a paired design if we have it and a
   weaker unpaired one if we do not.
3. **Pathway cliques.** Is a graph that is 90.6% MSigDB co-membership the graph
   we want? Options include capping clique size far below 100, weighting edges
   by set size, or dropping `protein_same_pathway` entirely and testing whether
   anything is lost.
4. **Correlation edge selection.** Should edge selection use a significance test
   with FDR control and a minimum overlapping-sample count, rather than a pooled
   80th-percentile cut on the raw correlation? Variant counts are already
   costed in `thesis/analysis/PLAN_improvement_experiment.md`.
5. **Scope of the re-preprocessing.** The DEPOD layer is 0.03% of edges and the
   map truncation costs 1.5% of PPI candidates. Neither changes a result. The
   pathway-clique share and the correlation-edge selection are the two choices
   that plausibly do. I would audit in that order rather than rebuild
   everything.

---

## 8. How to re-derive every number here

All read-only. None of these write to `src/`, `data/` or `outputs/`.

```bash
# edge composition of the generic graph
zcat outputs/edges.csv.gz | tail -n +2 | awk -F',' '{print $3}' | sort | uniq -c | sort -rn

# protein roles
zcat outputs/nodes.csv.gz | awk -F',' 'NR==1{for(i=1;i<=NF;i++) if($i=="protein_role") c=i; next} {n[$c]++} END {for (k in n) print n[k], k}'

# overlapping-sample distribution of the fold-change correlation edges
zcat outputs/Liver_network_edges.csv.gz | awk -F',' '$3 ~ /^site_corr_fc_/ {c[$6]++; t++} END {for (n in c) printf "n=%s %d (%.2f%%)\n", n, c[n], 100*c[n]/t}'

# raw row counts
wc -l data/raw/*.csv data/raw/*.tsv data/raw/*.txt
zcat data/raw/PTMcode2_associations_between_proteins.txt.gz | wc -l

# gene map coverage of the PPI file
tail -n +2 data/raw/gene_protein.csv | awk -F',' '{a[$1]=1; a[$2]=1} END {print "map covers", length(a)}'
```

The DEPOD trace and the liver column listing were produced with short
throwaway pandas scripts that read the two files and printed counts. Both files
are read with `pd.read_excel` exactly as the pipeline reads them.
