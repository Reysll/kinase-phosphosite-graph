# Research log

One entry per experiment, newest last. Each entry records the date, the change,
why it was made, what it produced and what was decided. Numbers here are
measured. Anything not yet measured says so.

Branch: `research/context-model`, from `main`. The goal of this branch is a
model that works, not a documented thesis. Earlier negative results are in
`STATE.md` (findings 1 to 37) and `EMBEDDING_DEFECT.md`.

---

## 2026-10-27 | Step 0: fast screening harness

**Change.** Added `thesis/analysis/screen.py`. It runs any model configuration
on a fixed subsample of the frozen leave-one-out trials and scores it against
five pass/fail criteria plus a kinase-centric metric. No project code was
modified.

**Why.** Every configuration measured so far cost between one and six hours, so
only four have ever been run. The binding constraint on finding a working model
is the number of configurations that can be rejected per evening, not the
precision of any single one.

### Design decisions and the measurements that forced them

**The subsample is exact, not approximate.** `run_leave_one_out` uses its
`positive_edges` argument for three separate things: the training matrix, the
site to true kinase map and the trial list. Subsampling that frame would shrink
the training set, which changes the model being screened. The harness instead
monkeypatches `leave_one_out._prepare_trials`, following the precedent in
`null_network.py` and `random_negatives.py`, so `positive_edges` stays the full
6,909 edges for training while only the trial list is cut. Each trial fits its
own model from the same full training matrix and the embedding is fit once
outside the loop, so every retained trial produces the rank it would in a full
run. Verified: 60 of 60 probe trials reproduced the full liver FC run's
`adjusted_held_out_rank` exactly. The same argument holds for the predict-all
pass, which receives `positive_edges` unrestricted exactly as the shipped
runner does, so the global model and every masked model are unchanged and
`target_sites` only decides which sites get scored.

This exactness is what makes the subsample cheap to validate. The current
model's scorecard can be computed by restricting the existing full runs, with
nothing to run.

**The trial subsample had to become two phase.** The frozen set is 85% rich
(5,872 rich, 836 average, 201 poor). A proportional 1,000-trial draw leaves 29
poor trials, which is too few for criterion 1, which tests poor and average
separately. The first design was therefore balanced: all 201 poor, 350 average,
500 rich. That design failed its own validation. The context effect lives almost
entirely in rich kinases, so a balanced draw diluted it from the full run's
rank-biserial 0.380 down to 0.072, and 3 of 12 scorecard verdicts disagreed with
the full run. The fix is two overlapping phases drawn in one pass and covered by
one run: a proportional phase of 1,000 trials that drives every overall test,
and a boost phase that tops the thin strata up to 201 poor and 350 average for
the per-stratum tests. Union 1,401 trials, flagged by `in_prop` in
`out/screen/screen_trials.csv`. After the fix, 12 of 12 verdicts agree with the
full run.

**The site subsample is two sets for one predict-all pass.** 300 FC-valid sites
drawn at random carry criterion 4, because eta-squared is a property of the site
population and must not be measured on a hand-picked set. A random 300-site draw
would leave almost no kinase with 5 known substrates, so the kinase-centric
metric gets its own set: every substrate site of the 29 kinases that have at
least 5 known substrates among the FC-valid sites, which is 337 sites. The union
is 612 sites, scored once. Seed 20261027, saved to `out/screen/screen_sites.csv`.

**Thresholds.** Criteria 1, 2, 3 and 5 are relative, so their thresholds follow
from the design: the model must win with Holm-adjusted p below 0.05, and for
criterion 5 the effect must exceed the measured node2vec seed-noise ceiling of
|r| = 0.373. Criterion 4's two numbers are judgement calls, set just beyond what
node2vec already reaches (eta-squared 0.990 to 0.992, median per-kinase rank SD
8.2 to 10.2), so passing means being more site-aware than any model measured so
far. They are constants at the top of `screen.py` and should be reviewed.

### Reference scorecard: the current model (original spectral)

Computed from the full runs on disk, restricted to the subsample. Full-run
verdicts in brackets where they are reported separately in
`out/screen/reference_scorecard.txt`.

| Criterion | Generic | Liver FC |
|---|---|---|
| 1 beats popularity (overall, poor, average) | FAIL | FAIL |
| 2 beats degree baseline | FAIL | FAIL |
| 3 beats its own degree-preserving null | N/A | PASS |
| 4 site dependence | FAIL | FAIL |
| 5 context effect exceeds seed noise | PASS | PASS |
| 6 kinase-centric AUC above chance (added) | FAIL | PASS |

Criterion 1 fails overall and on average kinases but passes on poor kinases in
both networks (rank-biserial +0.626 generic, +0.480 liver FC), which reproduces
finding 36. Criterion 3 is N/A for Generic because the null randomizes the
correlation edges and Generic has none. Criterion 4 fails by a wide margin,
eta-squared 0.9998 and median per-kinase rank SD 1.02 on the 300 random sites.
Generic covers only 29 of those 300 sites, which is the known 12% coverage of
the small generic graph file and not a new result.

### The kinase-centric metric is the one place the model shows signal

Ranking sites for each kinase, rather than kinases for each site, is the
direction in which the popularity baseline is constant and therefore cannot
compete. Measured on the 337-site pool, 29 kinases with at least 5 substrates
each:

| Network | Median AUC | Mean AUC | Above 0.5 | Median substrate rank |
|---|---|---|---|---|
| Generic | 0.453 | 0.480 | 41% | 64.2% of pool |
| Liver FC | 0.567 | 0.576 | 76% | 42.1% of pool |

Liver FC is above chance (Wilcoxon p = 0.011) and Generic is not (p = 0.80).
Paired across the same 29 kinases, Liver FC beats Generic on 23 of 29, median
AUC difference +0.063, Wilcoxon p = 0.0002. The per-kinase spread is large and
not degenerate: RAF1 0.951, PRKD1 0.875, PAK1 0.802, PRKAA2 0.795, PRKCD 0.775
at the top, against CDK2 0.305, SGK1 0.316, CDK1 0.360 at the bottom. Full table
in `out/screen/reference_kinase_centric_liver_fc.csv`.

Two observations, neither yet a claim. First, several of the high-AUC kinases
(PRKD1, PRKCD, PAK1, PRKG1, RPS6KA1) are among the 12 kinases that KSEA finds
significantly regulated and that never reach top-1, top-5 or top-10 in the
site-centric leave-one-out results. The kinase-centric direction may therefore
be recovering signal that hub bias hides in the site-centric direction. Second,
AUC trends downward with substrate count (Spearman -0.295 on liver FC) but not
significantly (p = 0.12), so the mirror-image-of-hub-bias reading is suggestive
only and is not established on 29 points.

### Runtime

Measured on liver FC, spectral, this machine. The leave-one-out cost is a
two-point fit at 40 and 120 trials, because dividing a single small probe by its
trial count amortizes the fixed setup over too few trials and overestimates a
larger run.

| Component | Cost |
|---|---|
| one-time setup (graph load, embedding fit, training matrix) | 16.3 s |
| marginal cost per trial | 0.551 s |
| one LOO screen, 1,401 trials | 13 min |
| criteria 1, 2, 3, 5 (two networks plus the null) | 41 min |
| criterion 4, `--mode eta` | 0.3 min |
| criterion 6, `--mode infer` | 46 min |

A full screened configuration is therefore about 42 min for criteria 1 to 5 and
another 46 min for criterion 6, against roughly 4 h for the equivalent full
runs. Criterion 6 should be a second stage, run only on configurations that have
already cleared the first five.

**Why criterion 6 stays expensive.** The predict-all pass costs 46 min for 612
sites against 47 min for all 3,228, because its cost is set by the number of
known true edges on the target sites, and the kinase-centric site set is
substrate sites by construction. 662 of the 736 masked own-edge fits land on the
612-site subsample. Subsampling sites does not help here.

**Why criterion 4 does not need that pass.** A masked model is only required for
a pair that was a positive label in training. On the 300-site eta set those are
47 of 124,500 pairs, 0.04%. Measured: eta-squared and the median per-kinase rank
SD are identical to six decimal places whether the true pairs are included or
dropped. `--mode eta` therefore scores criterion 4 from the global model alone
and reproduces the full run exactly (eta-squared 0.999784, median SD 1.0186,
both matching; 124,472 of 124,500 ranks identical) in 18 s rather than 46 min.

### A defect found in the harness, and one it exposes in earlier work

The first version of the predict-all screen passed `positive_edges`
unrestricted, following `thesis/analysis/spectral_fixed_runs.py`. That is wrong.
`src_prediction/runners/run_inference_all.py` restricts `positive_edges` to the
FC-valid site set before passing it (line 145), so the global model trains on
736 positives rather than all 14,596. The screen was therefore not reproducing
the full run: only 65,592 of 253,980 ranks matched, and the median per-kinase
rank SD came out 0.29 against the full run's 1.02 on the same sites. Fixed by
restricting to the same 3,228 FC-valid sites, which `--mode eta` then confirmed
reproduces the original to six decimal places.

**This implicates the spectral null-space validation runs of 2026-10-27.**
`spectral_fixed_runs.run_infer` passes the unrestricted 14,596 positives, so the
two predict-all runs behind finding 35 differ from the originals in TWO ways,
`drop_null_space=True` and a different training positive set, not the one flag
the comparison claimed. Finding 35 reported the median per-kinase rank SD
falling from 1.18 to 0.24 under the fix. The unrestricted positive set alone
moves the same statistic from 1.02 to 0.29 on the eta sites with no embedding
change at all, so most or all of that reported drop is likely the confound
rather than the fix. The leave-one-out half of that validation is unaffected,
because `run_leave_one_out` was passed the frozen trials in both runs. Finding
35's eta-squared and rank-SD numbers need rerunning before they are cited.
Findings 33, 34, 36 and 37 do not depend on the predict-all pass.

**Decision.** The harness stands as the screening tool for this branch. The
reference scorecard is the bar: a configuration is worth a full run only if it
moves criterion 4, which is the criterion the current model fails worst and the
one that the two repaired defects both made worse. The kinase-centric metric is
promoted from a diagnostic to a reported metric on every future screen, because
it is the only direction in which the current model is above chance and the only
one the popularity baseline cannot win. Nothing beyond step 0 was run, and no
configuration other than the current model has been screened. Finding 35 is
flagged as needing a rerun; that rerun is not part of step 0 and has not been
started.

---

## 2026-10-27 | Why the generic and liver graphs report different protein counts

**Change.** Added `thesis/analysis/graph_protein_audit.py`, a read-only audit of
the protein and site populations of both graphs, and of each node table against
its own edge list. Output in `thesis/analysis/out/graph_protein_audit.txt`. No
project code, data or results were modified.

**Why.** Dr. Ayati asked why the two graphs hold different numbers of proteins
when they should be roughly equal.

### Result

The declared counts are already nearly equal, and that is the problem.

| | generic | liver | difference |
|---|---|---|---|
| protein nodes | 3,122 | 3,195 | +73 |
| site nodes | 10,196 | 13,410 | +3,214 |
| kinases | 420 | 420 | 0 |

The generic graph is a strict subset of the liver graph. No generic node is
dropped, which is expected, because the liver builder starts from the generic
node table and appends to it.

The liver edge list refers to **1,133 protein nodes that the liver node table
does not declare**. Every one of them appears as the source of a `has_site`
edge. The generic graph has no such omission, and neither graph has a site
without a parent edge. The true liver protein population is therefore 4,328,
roughly 1,200 more than generic rather than roughly equal to it. The apparent
near-equality of the declared counts is an artifact of the incomplete node
table.

### Cause

`src/graph/liver_network.py` draws the two node kinds from two different
sheets. Protein nodes come from the ProteinExpression sheet alone (line 135).
Site nodes and their `has_site` edges come from the Phosphorylation sheet (line
167), using that sheet's own gene column as the parent. A gene detected only in
phospho data therefore receives sites and a `has_site` edge but never a protein
node.

The asymmetry is large because the two sheets are very different sizes:

| sheet | rows | genes | genes passing `min_identified=6` |
|---|---|---|---|
| ProteinExpression | 444 | 439 | 95 |
| Phosphorylation | 6,571 | 2,310 | 2,226 |

Of the 95 surviving ProteinExpression genes, 22 were already generic protein
nodes, which leaves exactly the 73 observed additions. Of the 2,226 surviving
Phosphorylation genes, 1,313 are absent from the liver node table. 1,132 of
those emit at least one `has_site` edge. The remaining 181 produce no site at
all, because their modification strings do not parse to a single site.

### Consequence, which is structural and not only cosmetic

Lines 595 to 621 derive `protein_ids` from the node table and pass it as
`allowed_protein_ids` and `existing_protein_ids` to the MSigDB and PPI edge
builders. The 1,133 omitted proteins are therefore excluded from all pathway
and interaction wiring. Measured in the liver edge list:

| | median degree | mean | max |
|---|---|---|---|
| omitted proteins (n = 1,133) | 1 | 1.6 | 21 |
| declared proteins (n = 3,195) | 631 | 670.3 | 2,663 |

The only relation that touches them is `has_site`, 1,859 edges in total. So
**1,859 of 13,410 liver sites (13.9%) hang off a degree-1 leaf parent** with no
pathway or interaction context, instead of being wired into the graph through
their protein.

This is visible in the connectivity. With base relations only the liver
embedding graph has 998 connected components against generic's 34, and the
largest holds 15,099 of 17,731 nodes. Adding the FC correlation edges pulls it
down to 122 components, which reproduces the figure already recorded in
`CATCH_UP.md`. The correlation edges are doing the work of holding the graph
together that the missing protein wiring should be doing.

### A second, smaller defect found in passing

One of the 1,133 is the literal node `PROTEIN:nan`, carrying 21 sites named
`SITE:nan-<residue><position>`, which between them hold about 26,900
correlation edges. It comes from 25 phospho rows whose Gene Symbol is missing.
The builder guards with `if not gene`, but the gene has already been through
`str(...)`, so a missing value arrives as the string `"nan"`, which is truthy
and passes the guard. Four of the 21 ids are each fed by two different rows, so
measurements from unrelated proteins are averaged together under one meaningless
identity.

Checked and found sound: duplicate site ids are handled correctly. Across the
filtered phospho sheet 874 of 3,637 distinct site ids (24%) are fed by more than
one row, which is the expected result of a site being measured in several
peptides, and `_prepare_site_fold_change_matrix` collapses them with
`groupby(level=0).mean()` at line 236.

### Effect on existing results

None of the reported numbers are wrong because of this. The node table is read
in only two places in the prediction layer, `evaluation_set.py` lines 20 and 34,
for liver site nodes and for generic kinase nodes. Both populations are
complete, 13,410 sites with no omission and 420 kinases in both graphs. The
embedding graph is built from the edge list by `build_graph_from_edges`, which
creates nodes implicitly, so the 1,133 proteins do carry embeddings and do
participate in every run. What is wrong is the graph structure they participate
with, and the protein count the node table reports.

### A number the planned rebuild will need

The new direction includes restricting candidate kinases to those present in
protein expression. On the current data that filter is not usable:

| candidate kinases | count |
|---|---|
| in the graph | 420 |
| in ProteinExpression, unfiltered | 20 |
| in ProteinExpression, passing `min_identified=6` | 2 |
| phospho-detected, passing `min_identified=6` | 106 |

The 20 are CAMK2A, CHEK2, ENPP3, ERBB2, GSK3A, LIMK2, MAP3K7, MAPK8, MASTL,
MYO3A, NEK6, PIK3C2A, PIK3CB, PLK1, RIPK2, SGK3, SMG1, TAOK2, TTK and ULK3.
Restricting to protein expression as it stands would cut the candidate set from
420 to 20, or to 2 under the current identification filter. This is consistent
with the note in `CLAUDE.md` that protein co-expression edges are disabled
because only 3 rows are usable. Either the filter has to mean phospho-detected,
which gives 106 of 420, or it depends on the larger protein expression table
that is not yet available.

**Decision.** Nothing fixed, because the fix is in `src/graph/liver_network.py`
and project code is not to be modified. Two changes are queued for the rebuild
that Dr. Ayati has asked for: register every `has_site` parent as a protein node
before `protein_ids` is handed to the MSigDB and PPI builders, and reject a
missing gene symbol rather than letting it through as the string `"nan"`. Both
change the graph, so both invalidate every result built on the current one, which
is acceptable only because the rebuild is already planned. The audit script is
read-only and can be rerun on the rebuilt graph as a regression check.
