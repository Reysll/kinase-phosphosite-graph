# Advisor brief — validation findings and a decision request

[NAME] · [DATE]

Every number comes from `thesis/analysis/out/`. No project code, data or results
were modified to produce them.

## What I tested

Before writing up the claim that liver context improves kinase-substrate
prediction, I ran five checks that had never been done: a significance test, a
random-network null, a seed-stability test, a test of whether the model uses the
site it is scoring, and a comparison against a trivial baseline.

## What each test showed

**1. Significance.** At trial level every context network beats Generic and
spectral beats node2vec, all significant after Holm correction (n = 6,894;
spectral median 79 → 74, r = +0.380). That part holds.

**2. Degree-preserving null.** Randomizing *which* sites the Liver FC
correlation edges connect, while preserving every site's edge count exactly,
reproduces the real network in 95.8–96.2% of trials and reproduces **100%** of
the necroptosis cluster's published rank shift.

**3. Seed stability (node2vec).** Seeds 7 and 13 move the Generic median from
188 to 177 to 148; the largest seed-only effect (+0.373) exceeds the context
effect (+0.249), seed 13 reverses the direction, and all six runs produced a
different dominant top-1 kinase.

**4. Site-dependence.** Kinase identity alone explains **99.97%** of rank
variance under spectral (η² = 0.9997). A typical kinase moves about **1 rank
position** across all 3,228 sites, and 48% move ≤ 1 rank anywhere. The model
applies a near-fixed popularity ordering rather than scoring the site.

**5. Popularity baseline — the decisive test.** I ranked candidates by substrate
count alone: no embedding, no training, no seed. Parity was verified first and
the baseline rebuilt to match the model exactly — the model trains only on the
6,909 frozen multi-kinase edges (not all 14,596 PSP edges) and scores 415
candidates (not 420). Correcting both made the baseline *slightly stronger*.

| Method | Median | Hit@10 | Hit@20 |
|---|---|---|---|
| **Popularity (parity-corrected)** | **25.0** | **33.4%** | **45.1%** |
| Liver FC spectral | 74.0 | — | 22.1% |
| Generic spectral | 79.0 | — | 10.6% |
| Liver FC / Generic node2vec (seed-avg) | 158.3 / 174.3 | — | — |

Holm-corrected, all p < 1e-300, popularity wins every aggregate comparison
(r = −0.686 to −0.873).

**6. Stratified — and here the models do win.**

| Stratum | Rank type | Winner | Effect |
|---|---|---|---|
| **Poor** (198 trials) | global | **models** | r +0.48 to +0.81, sig. |
| Poor | within-category | popularity | r −0.84 to −0.95 |
| Average (824) | global, trial | popularity | r −0.24 to −0.32 |
| **Average** (128 kinases) | global, kinase | **tie** | p = 1.00, n.s. |
| Average | within-category | popularity | r −0.85 to −0.95 |

For poorly annotated kinases the models place the true kinase far higher than
substrate count can (median 303 → 192–231). Within category, popularity still
orders better everywhere, which is what stratification was meant to fix.

**7. Liver-measured subset.** On the 510 trials whose site is in the liver
proteomics, the spectral context effect persists (+0.252), node2vec does not
reach significance, and the null-vs-real contrast stops being significant —
though at 92% ties this is underpowered, not evidence of equivalence.

## Three causes

1. **Site-independent popularity prior (confirmed).** η² = 0.9997, so
   correlation edges can only shift a global ordering. A high η² alone would not
   condemn the model, since real kinase activity does have a popularity
   component; the baseline is what settles it, and the model reproduces a
   *worse* prior than counting substrates. The exception is poor kinases, where
   the embedding does contribute.
2. **Low-sample correlation edges.** Only 8.0% of the 976,855 FC edges use the
   full n = 18; 92% use fewer, most at the floor of n = 6. (Mean |r| does not
   decline with n, so this is estimator variance, not inflation.)
3. **node2vec seed instability.** Seed noise exceeds the measured effect.

## Claims to retract or reframe

- **Poster/README node2vec numbers** (188 vs 165–174) — single seed; cross-seed
  spread is larger than the difference.
- **Hub-monopoly tables** (ULK1 100%, VRK2 94.8%, VRK3 95.5%) — describe seed
  42, not the method; six runs gave six different winners.
- **Necroptosis cluster as liver-specific** — the degree-matched random network
  reproduces the entire shift. The kinases are genuinely HCC-associated in the
  literature, but this analysis does not show the signal comes from liver
  co-variation.
- **"Spectral halves the median rank"** — true at trial level, but confined to
  well-annotated kinases and not significant per kinase.

**What survives:** graph construction and integration, the leakage fix, the
adjusted-rank metric, the four-network design, the spectral-vs-node2vec
trial-level comparison, KSEA (12/64 significant two-tailed), and one genuine
positive — for poorly annotated kinases the models beat substrate count on
global rank (median 303 → 192–231, r +0.48 to +0.81, all significant). That is
evidence the graph contributes real information for exactly the kinases the
project set out to help, even though the aggregate metric hides it.

One structural caveat: PSP-based leave-one-out is 85.2% rich-kinase trials
(5,872 of 6,894), so the headline metric structurally rewards annotation bias.
This mirrors the published imbalance in which ~80% of annotated phosphosites
belong to ~20% of kinases (Maier et al., 2025, *Mol Cell Proteomics* 24(8),
100994, https://doi.org/10.1016/j.mcpro.2025.100994).

## Decision I need

**A — report as an evaluation thesis.** Pipeline plus validation findings as the
contribution; fixes become future work. **0 h compute**, all results exist.

**B — run improvement experiments first.** Edge-set variants, each against its
own null (plan: `thesis/analysis/PLAN_improvement_experiment.md`): (a) FDR
significance, 846,586 edges; (b) top 5%/1% |r|, 244,214/48,843; (d) n ≥ 12 and
n ≥ 12 + top 5%, 379,783/77,356; (c) differential correlation, needs
recomputation. 12 runs, **~14–15 h**.

**C — middle path.** Variant (d) `n ≥ 12` only as a data-quality check,
**~2.1 h**, everything else as A.

**D — target the prior directly.** Pair features are `[K, S, |K−S|, K*S]`; the
raw `K` block lets the classifier learn a per-kinase intercept, i.e. a
popularity prior. Training on relational features only, `[|K−S|, K*S]`, attacks
the diagnosed cause rather than the edge set. **Not yet costed** — needs a
thesis-side variant of `build_pair_feature_table`, a decision on how many
networks to re-run, and its own null. Earlier exploratory work on this may exist
outside the repository; I have not seen it and report no result from it.

I recommend **A or C**, with **D** as the most scientifically direct next step
if there is appetite for compute. The judgement I would value most: are the
site-independence and popularity-baseline results a limitation of this pipeline,
or the thesis's central finding?
