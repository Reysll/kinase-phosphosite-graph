# CHAPTER V: EXPERIMENTS AND RESULTS

This chapter reports what the experiments measured. It does not argue for an
interpretation; Chapter VI does that.

Every quantity is traceable. Each table and claim carries an HTML comment naming
the output file it came from, and `thesis/STATE.md` lists the command that
regenerates each file. Effect sizes accompany every p-value, because several
results in this chapter combine very small p-values with effect sizes close to
zero, and reporting only one of the two would misrepresent them.

Unless stated otherwise, all leave-one-out metrics use **n = 6,894** trials, the
6,909 frozen trials minus the 15 that are unrankable in every run, and all
rankings cover the **415** candidate kinases that receive embeddings.

Two conventions hold throughout. **Lower rank is better.** For paired tests the
difference is `rank_A - rank_B`, so a **positive** median paired difference and a
**positive** matched-pairs rank-biserial correlation both mean configuration B
achieved the better rank.

## 5.1 The Four-Network Comparison

Four networks were compared under each embedding: Generic with no correlation
edges, and Control, Cancer and Liver FC with their respective correlation edge
sets. Ten paired Wilcoxon signed-rank tests form one family, Holm-corrected
together.

### 5.1.1 Trial-Level Results

<!-- source: thesis/analysis/out/stat_tests.txt, stat_tests_trial_level.csv -->

| Embedding | Comparison | Median A | Median B | Rank-biserial | p (Holm) | Favours |
|---|---|---|---|---|---|---|
| node2vec | Generic vs Control | 188.0 | 166.5 | +0.175 | 5.53e-36 | Control |
| node2vec | Generic vs Cancer | 188.0 | 174.0 | +0.037 | 0.0085 | Cancer |
| node2vec | Generic vs Liver FC | 188.0 | 165.0 | +0.249 | 3.21e-70 | Liver FC |
| spectral | Generic vs Control | 79.0 | 74.0 | +0.391 | 1.95e-172 | Control |
| spectral | Generic vs Cancer | 79.0 | 74.0 | +0.380 | 6.77e-163 | Cancer |
| spectral | Generic vs Liver FC | 79.0 | 74.0 | +0.380 | 1.40e-162 | Liver FC |

All six context comparisons are significant after Holm correction. Every
context network outranks Generic at the trial level under both embeddings.

The effect sizes are not uniform. Under node2vec the Cancer comparison has a
rank-biserial of +0.037, which is significant only because n is 6,894. Under
spectral the three context networks are nearly indistinguishable from each
other, all reaching a median of 74.0.

### 5.1.2 Kinase-Level Results

Aggregating to one value per held-out kinase, taking the mean adjusted rank per
kinase and testing over the 355 kinases that appear as held-out kinases, gives a
different picture.

<!-- source: thesis/analysis/out/stat_tests_kinase_level.csv -->

| Embedding | Comparison | Median A | Median B | Rank-biserial | p (Holm) | Significant |
|---|---|---|---|---|---|---|
| node2vec | Generic vs Control | 198.4 | 200.9 | +0.022 | 1.000 | no |
| node2vec | Generic vs Cancer | 198.4 | 202.5 | +0.008 | 1.000 | no |
| node2vec | Generic vs Liver FC | 198.4 | 201.0 | +0.007 | 1.000 | no |
| spectral | Generic vs Control | 197.0 | 196.5 | **-0.232** | 0.0012 | yes |
| spectral | Generic vs Cancer | 197.0 | 197.5 | **-0.247** | 5.66e-04 | yes |
| spectral | Generic vs Liver FC | 197.0 | 197.0 | **-0.242** | 8.12e-04 | yes |

Under node2vec no context comparison survives correction. Under spectral all
three are significant and the sign is **negative**, meaning they favour
**Generic**.

These two sections answer different questions and both answers stand. The
trial-level test asks whether a randomly chosen held-out edge is ranked better
under context, and the answer is yes. The kinase-level test asks whether a
randomly chosen kinase is ranked better under context, and under spectral the
answer is no. The difference arises because trials are not evenly distributed
across kinases: 5,872 of the 6,894 trials involve rich kinases, so the
trial-level result is dominated by them.

### 5.1.3 How Similar Are the Context Networks?

<!-- source: thesis/analysis/out/t12_spearman.csv, t12_top10_overlap.csv -->

| Embedding | Pair | Spearman rho | Fraction of trials with identical rank |
|---|---|---|---|
| spectral | Control vs Cancer | 0.99937 | 56.7% |
| spectral | Control vs Liver FC | 0.99953 | 56.5% |
| spectral | Cancer vs Liver FC | 0.99929 | 65.2% |
| spectral | Generic vs Control | 0.8679 | 0.6% |
| node2vec | Control vs Liver FC | 0.7135 | 0.9% |

Under spectral the three context networks are almost the same ranking. Their
per-site top-10 lists share a mean of **9.57 of 10** kinases for Control against
Liver FC, and **57.0%** of the 3,228 sites have an identical top-10 under
Control and Liver FC. Which correlation type is used barely changes the output.
Generic, by contrast, is clearly distinct from all three.

**Proposed Table 5.1**: the trial-level and kinase-level comparisons side by
side, as above.
**Proposed Figure 5.1**: `outputs_prediction/analysis/poster_loo_figure.png`,
median adjusted rank by network and embedding.

## 5.2 The Embedding Comparison

The same family of tests contains four embedding comparisons, holding the
network fixed.

<!-- source: thesis/analysis/out/stat_tests.txt, block B -->

| Network | node2vec median | spectral median | Rank-biserial | p (Holm) |
|---|---|---|---|---|
| Generic | 188.0 | 79.0 | +0.508 | 1.08e-291 |
| Control | 166.5 | 74.0 | +0.626 | 0.000 |
| Cancer | 174.0 | 74.0 | +0.640 | 0.000 |
| Liver FC | 165.0 | 74.0 | +0.583 | 0.000 |

At the trial level spectral roughly halves the median rank in every network,
with large effect sizes.

At the kinase level, **none of the four embedding comparisons is significant**
(all p Holm = 1.000, rank-biserial +0.015 to +0.049).

<!-- source: thesis/analysis/out/stat_tests_kinase_level.csv -->

So the embedding advantage has the same structure as the context advantage: it
is strong per trial and absent per kinase. The qualification matters for how the
result is stated. "Spectral halves the median rank" is accurate at the trial
level and not demonstrated per kinase.

Restricting to kinases whose rank actually varies across sites does not change
the conclusion. Under that restriction the trial-level embedding effects remain
significant with rank-biserial +0.354 to +0.545.

<!-- source: thesis/analysis/out/stat_tests.txt, section 3; stat_tests_informative.csv -->

**Proposed Table 5.2**: as above, trial level and kinase level together.

## 5.3 Seed Variance

This analysis applies to node2vec only, for the reason given in Section 4.9.2:
spectral is deterministic given its `random_state`, so repeating it returns
identical output. Three seeds were run, 42, 7 and 13.

<!-- source: thesis/analysis/out/seed_averaged.txt, seed_averaged_medians.csv -->

| Network | Seed 42 | Seed 7 | Seed 13 | Spread | SD of medians |
|---|---|---|---|---|---|
| Generic | 188.0 | 177.0 | 148.0 | **40.0** | 20.66 |
| Liver FC | 165.0 | 157.0 | 153.0 | 12.0 | 6.11 |

The Generic median moves 40 rank positions across seeds. Per trial, the standard
deviation of a trial's rank across the three seeds has a mean of **42.4** in
Generic and **39.1** in Liver FC, with maxima of 186.8 and 196.2.

The decisive comparison places the seed noise next to the context effect
measured in Section 5.1.

<!-- source: thesis/analysis/out/seedvar_compare.txt, seedvar_compare.csv -->

| Contrast | Type | Median paired diff | Rank-biserial |
|---|---|---|---|
| Generic vs Liver FC, seed 42 | context effect | 7.0 | **+0.249** |
| Generic seed 42 vs seed 13 | seed noise | 17.0 | **+0.373** |
| Generic seed 42 vs seed 7 | seed noise | 4.0 | +0.071 |
| Liver FC seed 42 vs seed 13 | seed noise | -3.0 | -0.020 |

Changing the random seed alone produces a larger effect (+0.373) than changing
the network (+0.249), and a larger median shift (17.0 against 7.0). One
seed-to-seed contrast, Generic 42 against 13, also reverses relative to the
direction the context comparison reports.

Hub identity is not stable either.

<!-- source: thesis/analysis/out/seed_averaged_top1.csv -->

| Network | Seed 42 | Seed 7 | Seed 13 |
|---|---|---|---|
| Generic | ULK1 (76.5%) | YES1 (98.7%) | WNK4 (80.9%) |
| Liver FC | ULK1 (100%) | WEE1 (96.7%) | ZAP70 (91.2%) |

Six runs produced six different most-frequent top-1 kinases, each dominating a
large share of trials within its own run.

**Proposed Table 5.3**: per-seed medians and the noise-versus-effect comparison.
**Proposed Figure 5.2**: new, a dot plot of per-seed median rank per network with
the context effect drawn as a bracket, making the noise-to-effect ratio visible.

## 5.4 The Degree-Preserving Null

The null replaces the Liver FC correlation edges with random edges of identical
degree sequence, destroying about 99.9% of the original pairings.

### 5.4.1 Leave-One-Out

<!-- source: thesis/analysis/out/null_compare.txt, null_compare.csv -->

| Comparison | Median A | Median B | Rank-biserial | Ties | p (Holm) |
|---|---|---|---|---|---|
| Generic vs Null (seed 1) | 79.0 | 74.0 | +0.380 | 0.70% | 7.05e-163 |
| Generic vs Liver FC | 79.0 | 74.0 | +0.380 | 0.70% | 7.05e-163 |
| Null (seed 1) vs Liver FC | 74.0 | 74.0 | +0.264 | **95.75%** | 9.22e-06 |

The first two rows are identical to three decimal places. The randomized network
improves on Generic exactly as much as the real one does.

Comparing the real and null networks directly, **95.75% of trials receive an
identical rank** under seed 1 and **96.17%** under seed 2. The remaining
difference is significant (p Holm = 9.22e-06, rank-biserial +0.264) but it is
computed over the 4.25% of trials that differ at all.

The agreement holds within every category:

| Category | n | Generic | Null | Liver FC |
|---|---|---|---|---|
| rich | 5,872 | 69.0 | **53.0** | **53.0** |
| average | 824 | 209.0 | **212.0** | **212.0** |
| poor | 198 | 231.0 | **272.5** | **272.5** |

Null and Liver FC are identical in all three. Both also make poor and average
kinases worse than Generic.

### 5.4.2 The Necroptosis Cluster

A cluster of IRAK4, RIPK1, RIPK3, MLKL and TRPM7 had been identified as gaining
predicted relevance under the disease networks, and RIPK1, RIPK3 and MLKL are
the core necroptosis kinases.

<!-- source: thesis/analysis/out/null_casestudy.txt, null_casestudy_medians.csv -->

| Kinase | Control | Liver FC | Null | Real shift | Null shift | Null reproduces |
|---|---|---|---|---|---|---|
| IRAK4 | 308.0 | 203.0 | 203.0 | -105.0 | -105.0 | **100%** |
| RIPK1 | 175.0 | 86.0 | 86.0 | -89.0 | -89.0 | **100%** |
| TRPM7 | 367.0 | 286.0 | 286.0 | -81.0 | -81.0 | **100%** |
| RIPK3 | 280.0 | 254.0 | 254.0 | -26.0 | -26.0 | **100%** |
| MLKL | 340.0 | 318.0 | 318.0 | -22.0 | -22.0 | **100%** |
| TTK | 6.0 | 7.0 | 6.0 | +1.0 | 0.0 | 0% |

The randomized network reproduces the shift of all five cluster kinases exactly,
to the rank. TTK is listed for completeness and is not evidence of anything: its
shift is a single rank position, which is too small to evaluate against the null.

### 5.4.3 Within-Category Ranking Against the Null

Kinases were re-ranked inside their own category at each of the 3,228 sites and
the real network compared with the null per kinase by KS test with BH-FDR.

<!-- source: thesis/analysis/out/withincat_null.txt, withincat_null_summary.csv -->

| Category | Tested | Differ from null (BH-FDR < 0.05) |
|---|---|---|
| poor | 162 | 48 (29.6%) |
| average | 132 | 62 (47.0%) |
| rich | 121 | 74 (61.2%) |
| **all** | **415** | **184 (44.3%)** |

Read alone, 44.3% of kinases differ from the null. The effect sizes contradict
that reading.

<!-- source: thesis/analysis/out/withincat_null.txt, effect-size section -->

| Quantity | Value |
|---|---|
| Of the 184 significant kinases, median shift exactly 0 ranks | **159** |
| Shift of 1 rank | 23 |
| Shift of 2 ranks | 1 |
| Largest shift anywhere | **2 ranks** |
| Pearson r, real shift against null shift, all 415 kinases | **0.9944** |
| Spearman, same | 0.9600 |
| Median gap between the two shifts | 0 |

Each kinase is tested across 3,228 sites, so the KS test resolves distributional
differences far below one rank position. For reference, Control against Liver FC
is significant for 75.9% of kinases and Control against Null for 74.9%, which
are the same figure within noise.

**Proposed Table 5.4**: the three null comparisons and the category breakdown.
**Proposed Table 5.5**: the necroptosis cluster as above.
**Proposed Figure 5.3**: `outputs_prediction/analysis/spotlight_necroptosis.png`,
with the null overlaid; the existing figure shows only real against Control and
needs the null series added.
**Proposed Figure 5.4**: new, a scatter of real shift against null shift across
all 415 kinases with the identity line, which makes r = 0.9944 immediate.

## 5.5 Site Dependence

<!-- source: thesis/analysis/out/site_dependence.txt, site_dependence.csv -->

| Embedding | Network | eta-squared | Median per-kinase rank SD | Kinases with SD <= 1 |
|---|---|---|---|---|
| spectral | Control | **0.9997** | 1.05 | 48.2% |
| spectral | Liver FC | **0.9997** | 1.18 | 42.2% |
| spectral | Generic | 0.9886 | 4.93 | 6.0% |
| node2vec | Control | 0.9903 | 9.08 | 3.4% |
| node2vec | Liver FC | 0.9920 | 8.20 | 3.9% |
| node2vec | Generic | 0.9879 | 10.20 | 2.4% |

Under spectral, kinase identity accounts for **99.97%** of the variance in
predicted rank across the 3,228 sites. The median kinase moves about **one rank
position** between sites, and 96% of kinases move 5 positions or fewer.

The scorer therefore applies a near-fixed global ordering of kinases and
conditions very weakly on which site it is scoring. node2vec is less extreme,
with a median per-kinase spread of 8 to 10 positions, but is still above 0.98.

**Proposed Table 5.6**: as above.
**Proposed Figure 5.5**: new, the distribution of per-kinase rank standard
deviation for spectral and node2vec, with the 1-position mark highlighted.

## 5.6 The Popularity Baseline

The baseline ranks kinases by their substrate count, identically at every site,
using the same 6,909 frozen training edges and the same 415 candidates.

<!-- source: thesis/analysis/out/popularity_baseline.txt, popularity_parity_global.csv -->

| Method | Median adjusted rank | Hit@10 | Hit@20 |
|---|---|---|---|
| **Popularity baseline** | **25.0** | **33.42%** | **45.08%** |
| Liver FC spectral | 74.0 | | 22.1% |
| Generic spectral | 79.0 | | 10.6% |
| Liver FC node2vec, seed-averaged | 158.3 | | |
| Generic node2vec, seed-averaged | 174.3 | | |

<!-- source: thesis/analysis/out/popularity_parity_tests_all.csv -->

| Comparison | Median paired diff | Rank-biserial | p (Holm) | Favours |
|---|---|---|---|---|
| Popularity vs Generic spectral | -47.0 | **-0.758** | 0.000 | Popularity |
| Popularity vs Liver FC spectral | -26.0 | **-0.686** | 0.000 | Popularity |
| Popularity vs Generic n2v seed-avg | -126.3 | -0.873 | 0.000 | Popularity |
| Popularity vs Liver FC n2v seed-avg | -108.3 | -0.863 | 0.000 | Popularity |

A model-free baseline outperforms every trained configuration on the global
metric, with large effect sizes. Correcting an earlier parity error made the
baseline stronger, moving its median from 26.0 to 25.0.

### 5.6.1 Stratified Results

The aggregate hides a reversal by category.

<!-- source: thesis/analysis/out/popularity_parity_stratified.csv -->

Baseline performance by category, global rank:

| Category | Baseline median | Hit@20 |
|---|---|---|
| rich | 18.0 | 52.93% |
| average | 162.0 | **0.00%** |
| poor | 303.0 | **0.00%** |

For poor kinases, global rank, trial level:

| Comparison | Popularity | Model | Rank-biserial | p (Holm) | Favours |
|---|---|---|---|---|---|
| vs Generic spectral | 303.0 | 231.0 | **+0.626** | 5.34e-14 | model |
| vs Liver FC spectral | 303.0 | 272.5 | **+0.480** | 5.34e-09 | model |
| vs Generic n2v seed-avg | 303.0 | 194.7 | **+0.808** | 3.37e-22 | model |
| vs Liver FC n2v seed-avg | 303.0 | 192.2 | **+0.770** | 1.77e-20 | model |

Every model beats the baseline for poor kinases on global rank, with large
positive effect sizes, and the same holds at kinase level (+0.625 to +0.831).

The direction flips again under within-category rank, where the baseline wins
for poor kinases (rank-biserial -0.841 to -0.954, all p Holm < 1e-23), because
within-category ranking removes the cross-category competition that the baseline
handles badly.

**Proposed Table 5.7**: global comparison and the poor-kinase stratum together.
**Proposed Figure 5.6**: new, median rank by category for baseline and the two
spectral networks, grouped bars, making the aggregate-versus-stratified
divergence visible in one view.

## 5.7 KSEA and Its Relation to Ranking

<!-- source: thesis/analysis/out/ksea_twotailed.txt, ksea_twotailed.csv -->

KSEA was recomputed with the two-tailed p-value appropriate to a signed score,
`2 * norm.sf(|score|)`, in place of the one-tailed form used previously.

| Convention | Significant at BH-FDR < 0.05 |
|---|---|
| One-tailed, as previously published | 17 of 64 |
| **Two-tailed, used here** | **12 of 64** |

Five kinases do not survive the correction: CHEK1, MAPK3, PRKCE, PRKDC and TBK1.

The overlap between the KSEA set and the kinases the model ranks highest was
computed over the 64 kinases both methods can score.

<!-- source: thesis/analysis/out/ksea_ranking_overlap.txt, ksea_ranking_overlap.csv -->

| Ranking cutoff | Ranking set size | Overlap | Expected by chance | Depletion p |
|---|---|---|---|---|
| top 10 | 3 | 0 | 0.56 | 0.530 |
| top 20 | 8 | 1 | 1.50 | 0.533 |
| top 50 | 19 | 3 | 3.56 | 0.495 |

The observed overlap is essentially the expected overlap at every cutoff, and no
hypergeometric test of depletion approaches significance in any of the six
comparisons across the two context networks. The two methods nominate few
kinases in common, but no more rarely than chance. The small overlap follows
from both nominating small sets from a pool of 64, not from systematic
disagreement.

**Proposed Table 5.8**: the two-tailed KSEA result and the overlap table.
**Proposed Figure 5.7**: `outputs_prediction/analysis/ksea_vs_loo_rank.png`,
KSEA score against LOO rank, updated to mark two-tailed significance.

## 5.8 Differential Phosphosite Analysis

This analysis uses the liver phosphoproteomics directly and involves neither the
graph nor the ranking model.

### 5.8.1 Normalization

<!-- source: thesis/analysis/out/diffsites/normalization_check.txt, channel_stats.csv -->

The `Phosphorylation` sheet is labelled "Abundance" rather than "Abundances
(Normalized)", and the measurement confirms it is not normalized. Per-channel
medians span a factor of **6.13**, or **2.62 log2 units**, from 16.60 to 19.22.
Per-channel median centering on log2 was applied before any test, reducing the
spread to 0.000.

The channel medians also reveal the batch structure. Within a batch the three
same-condition channels agree to about 0.15 log2, while across batches they
differ by up to 2.62, and the tumour-against-control offset changes direction
between batches.

### 5.8.2 Testable Set

<!-- source: thesis/analysis/out/diffsites/diffsites.txt -->

| Stage | Count |
|---|---|
| Rows in the sheet | 6,571 |
| Parsing to exactly one phosphosite | 5,249 |
| Unique sites after collapsing peptides | 3,780 |
| **Tested, >= 6 observed in both groups** | **2,973** |
| Strict tier, >= 12 observed in both groups | 2,046 |

Overall missingness is 34.7%. No imputation was applied.

### 5.8.3 Primary and Sensitivity Tests

| Test | Significant at BH-FDR < 0.05 |
|---|---|
| **Linear model, condition + batch (primary)** | **1,548 of 2,973 (52.1%)** |
| Welch unpaired, batch ignored | 939 (31.6%) |
| Mann-Whitney U | 861 (29.0%) |

Of the 1,548 significant under the primary test, 801 are up in tumour and 747
down. 783 have |log2FC| >= 1, and **659** meet the shortlist criteria of
BH-FDR < 0.05, |log2FC| >= 1 and at least 9 observed values per group.

The batch term matters:

| Outcome | Sites |
|---|---|
| Significant under both tests | 922 |
| **Gained by modelling batch** | **626** |
| Lost by modelling batch | 17 |
| **Changing significance in total** | **643 (21.6% of tested)** |

The gain comes from reduced residual variance rather than shifted estimates:
raw and batch-adjusted log2 fold changes correlate at **r = 0.9649**. The strict
tier behaves consistently, with 1,128 of 2,046 sites (55.1%) significant.

The strongest changes are coherent with liver biology, the largest being
CYP2A6-S131 at -5.52 log2 and PGM1-S509 at -5.36, with SQSTM1 and SPP1 sites
among the largest increases.

### 5.8.4 Presence and Absence

114 sites are detected in one group and never in the other, 77 tumour-only and
37 control-only, and all 114 reach significance by Fisher exact test on
detection counts. As Section 5.8.6 notes, that figure is largely a property of
the selection filter.

<!-- source: thesis/analysis/out/diffsites/diffsites_limitations.txt -->

Their batch distribution was checked separately: **none** of the 114 has its
detections confined to a single batch, 66.7% span four or more of the six
batches, and for **all 114** every batch containing a detection also had the
opposite group well covered. The one-sided detection is therefore not an
artifact of a batch in which the other condition was simply not measured.

### 5.8.5 Kinase Mapping

<!-- source: thesis/analysis/out/diffsites/diffsites.txt, diffsites_kinase_map.csv, diffsites_ksea.csv -->

Direct PhosphoSitePlus annotation covers a small minority of sites: **372 of
2,973 tested sites (12.5%)** have at least one annotated kinase, and **209 of
1,548 significant sites (13.5%)**. Unannotated sites are retained in every
output rather than filtered out.

KSEA run on the differential result, two-tailed, scores 64 kinases with at least
three measured substrates and finds **15 significant**. The largest scores are
PRKACA at -5.68 and RPS6KA1 at -3.66 downward, CDK2 at +4.58 and CSNK2A1 at
+3.32 upward.

All 15 KSEA-significant kinases also appear among the 109 kinases annotated to
at least one significant site. That containment is close to structural and is
not evidence of agreement between independent methods, for the reason given in
Section 5.8.6.

### 5.8.6 Limitations of This Analysis

<!-- source: thesis/analysis/out/diffsites/diffsites_limitations.txt -->

Four limitations qualify everything in Section 5.8.

**No correction for protein abundance.** A phosphosite can change because the
protein carrying it changed. The `ProteinExpression` sheet cannot support a
correction: it holds 444 rows, is 91.7% empty, has a median of 0 observed values
per row, and only **2 rows** have at least 6 observations in both groups. Only
**81 of 2,973 tested sites (2.7%)** and **19 of 659 shortlist sites (2.9%)**
have their protein present in that sheet at all, and **0 of 659** have protein
data sufficient for a correction. Its caption also describes it as proteins
unique to tumour tissue rather than a matched quantitative matrix.

**Presence and absence sites are candidates, not a statistical finding.** They
are selected by a detection-count filter and then tested with a detection-count
test, so all 114 reaching significance reflects the filter.

**The containment of the KSEA set within the annotation set is near-automatic.**
Route A collects every kinase annotated to any of 1,548 significant sites and
contains 109 kinases; KSEA scores only kinases with at least three measured
substrates, largely the same well-annotated kinases.

**One cohort, unpaired, unresolved provenance.** There is no external
replication, the pairing is not recorded so the test is unpaired, which is
conservative rather than anti-conservative, and the originating study is still
unidentified.

**Proposed Figure 5.8**:
`thesis/analysis/out/diffsites/diffsites_volcano.png`, batch-adjusted log2 fold
change against -log10 BH-FDR with the shortlist marked.
**Proposed Figure 5.9**: new, per-channel log2 medians before and after
centering, grouped by batch and condition, which shows both the normalization
problem and the batch structure in one panel.
**Proposed Table 5.9**: the three tests and the significance-switching counts.
**Proposed Table 5.10**: the top 30 shortlist sites with effect sizes, observed
counts and annotated kinases.

## 5.9 Negative Sampling

An audit of the training code established that negative examples are not sampled
but taken as a fixed prefix of an alphabetically sorted candidate list.

<!-- source: thesis/analysis/out/negative_sampling_audit.txt -->

| Quantity | Value |
|---|---|
| Distinct kinases ever used as a negative | **55 of 420** |
| Candidate-list positions occupied | indices 0 to 54 |
| Negative pairs generated | 122,950 |
| Times a typical such kinase is used | 2,458 or 2,459 of 2,459 sites |
| Scored kinases never used as a negative | **361 of 415 (87.0%)** |

Whether this advantaged the never-negative kinases was tested directly.

<!-- source: thesis/analysis/out/negative_sampling_audit.csv -->

| Source | Metric | Ever-negative | Never-negative | p | Rank-biserial |
|---|---|---|---|---|---|
| Predict-all, Control | median rank | 213.5 | 207.0 | 0.740 | -0.028 |
| Predict-all, Liver FC | median rank | 215.5 | 207.0 | 0.698 | -0.033 |
| LOO Generic | adjusted rank | 240.0 | 191.5 | 0.060 | -0.170 |
| LOO Generic | within-category | 77.0 | 65.5 | 0.085 | -0.156 |
| LOO Control | adjusted rank | 206.0 | 196.8 | 0.569 | -0.052 |
| LOO Liver FC | adjusted rank | 206.0 | 197.3 | 0.563 | -0.052 |

**No comparison reaches significance at 0.05.** The largest gap is in the
Generic network, 48.5 rank positions favouring never-negative kinases, at
p = 0.060 across eight uncorrected comparisons.

The defect in the training design is established. A consequent bias in the
reported rankings is **not** established, and no claim of one is made. No
experiment was rerun with a corrected sampler, since that would break
comparability with every existing run against the frozen trial set.

**Proposed Table 5.11**: as above.

## 5.10 Summary of Results

| # | Result | Section |
|---|---|---|
| 1 | All six context comparisons significant at trial level; effect sizes +0.037 to +0.391 | 5.1.1 |
| 2 | At kinase level, node2vec context effects vanish and spectral ones favour Generic | 5.1.2 |
| 3 | The three context networks are nearly identical under spectral, rho >= 0.9993 | 5.1.3 |
| 4 | Spectral halves the trial-level median rank; not significant per kinase | 5.2 |
| 5 | node2vec seed noise (+0.373) exceeds the context effect (+0.249); six runs, six hubs | 5.3 |
| 6 | The degree-preserving null matches Liver FC in 95.8 to 96.2% of trials | 5.4.1 |
| 7 | The null reproduces 100% of the necroptosis cluster shift | 5.4.2 |
| 8 | Within-category ranking does not separate real from null; r = 0.9944 | 5.4.3 |
| 9 | Kinase identity explains 99.97% of rank variance under spectral | 5.5 |
| 10 | A substrate-count baseline beats every model globally (25.0 against 74.0) | 5.6 |
| 11 | Every model beats that baseline for poor kinases on global rank | 5.6.1 |
| 12 | Two-tailed KSEA gives 12 of 64, not 17 | 5.7 |
| 13 | KSEA and ranking overlap no less than chance predicts | 5.7 |
| 14 | 1,548 of 2,973 phosphosites differ between tumour and adjacent tissue | 5.8.3 |
| 15 | Modelling batch changes significance for 643 sites (21.6%) | 5.8.3 |
| 16 | Negatives come from 55 of 420 kinases; 87% of scored kinases never appear | 5.9 |
| 17 | No measurable rank consequence of that defect | 5.9 |

Results 1 and 2, and 10 and 11, are pairs that point in opposite directions at
different levels of aggregation. Both members of each pair are reported, and
neither is presented as superseding the other. Chapter VI takes up what they
mean together.
