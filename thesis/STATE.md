# Thesis state — recovery point

Last updated: 2026-10-26 (fourth pass). Advisor decisions recorded below;
see "Current plan" for the three tasks now in flight.
All paths relative to the repo root. All analysis scripts are read-only with
respect to `src/`, `src_prediction/`, `data/` and `outputs/`, and write only to
`thesis/analysis/out/`. Interpreter: `.venv312/Scripts/python`.

## Findings, with source files and commands

| # | Finding | Output file | Command |
|---|---|---|---|
| 1 | Trial-level: all 10 context/embedding contrasts significant (Holm, n=6,894). Spectral 79→74, r=+0.380 | `out/stat_tests.txt`, `stat_tests_trial_level.csv` | `python -m thesis.analysis.stat_tests` |
| 2 | Kinase-level: node2vec context n.s.; spectral context favours Generic (r −0.23 to −0.25) | `out/stat_tests_kinase_level.csv` | same |
| 3 | KSEA two-tailed: **12/64** significant (published 17/64 used `norm.sf(abs(z))`) | `out/ksea_twotailed.txt` | `python -m thesis.analysis.ksea_twotailed` |
| 4 | Integrity: identical trials/kinases/NaN across all 8 runs; 15 unrankable; 355 kinases | `out/verify_tests.txt` | `python -m thesis.analysis.verify_tests` |
| 5 | Stratified: rich gains, poor/average worsen under context | `out/verify_stratified_{rich,poor}_{trial,kinase}.csv` | same |
| 6 | Within-category: rich gain disappears (33→32 n.s.); poor/average significantly worse | `out/t11_*.csv`, `out/tier1_analyses.txt` | `python -m thesis.analysis.tier1_analyses` |
| 7 | Context networks near-identical under spectral (ρ=0.999; 57–99% identical top-10) | `out/t12_spearman.csv`, `t12_top10_overlap.csv` | same |
| 8 | Edge accounting: base 1,133,920; Liver FC +976,855 = 2,110,775. 77–79% of corr edges touch no annotated site | `out/t13_*.csv` | same |
| 9 | **Null reproduces Liver FC**: 95.8–96.2% identical ranks, 3 seeds | `out/null_compare.txt`, `null_compare.csv` | `python -m thesis.analysis.null_compare` |
| 10 | **Null reproduces 100% of necroptosis shift**. TTK was reported as a 0% exception; **retracted**, its Control to Liver FC shift was 1 rank (6 to 7), too small to evaluate against the null | `out/null_casestudy.txt`, `null_casestudy_ks.csv` | `python -m thesis.analysis.null_casestudy` |
| 11 | Seed noise (+0.373) exceeds context effect (+0.249); seed 13 reverses direction | `out/seedvar_compare.txt` | `python -m thesis.analysis.seedvar_compare` |
| 12 | Six runs, six different top-1 hubs; Generic median spread 40 ranks | `out/seed_averaged.txt`, `seed_averaged_top1.csv` | `python -m thesis.analysis.seed_averaged` |
| 13 | **η² = 0.9997** (spectral); median per-kinase rank SD ≈ 1 across 3,228 sites | `out/site_dependence.txt`, `site_dependence.csv` | `python -m thesis.analysis.site_dependence` |
| 14 | **Popularity baseline beats all models** (median 25.0 vs 74.0; Hit@20 45.1% vs 22.1%) | `out/popularity_baseline.txt` | `python -m thesis.analysis.popularity_baseline` |
| 15 | **Models beat popularity for poor kinases, global rank** (r +0.48 to +0.81) | `out/popularity_parity_stratified.csv` | same |
| 16 | Liver-measured subset (n=510): spectral context survives; null contrast n.s. | `out/subset_analysis.txt` | `python -m thesis.analysis.subset_analysis` |
| 17 | **Within-category ranking does not escape the null.** 184/415 kinases (44.3%) differ from the null by KS+BH, but 159 of those 184 have a median shift of exactly 0 and the max is 2 ranks; liverfc_shift vs null_shift r=0.9944 | `out/withincat_null.txt`, `withincat_null_*.csv` | `python -m thesis.analysis.withincat_null` |
| 18 | **Identifier map truncated at 2^20 lines.** `gene_protein.csv` has exactly 1,048,576 lines (Excel's row cap); 1,664/16,185 PPI proteins unmapped (10.3%), costing 6,960/472,000 PPI rows (1.5%) | `GRAPH_WALKTHROUGH.md` §3.5 | see that file §8 |
| 19 | Liver sheets are supplementary Tables S4 and S11 of an unidentified HCC paper; 18+18 is consistent with 6 batches of 3+3 but the patient pairing is not recorded | `GRAPH_WALKTHROUGH.md` §5.1-5.2 | same |
| 20 | **Phospho sheet is NOT normalized**: per-channel medians span 6.13x (2.62 log2). Per-channel median centering applied before testing | `out/diffsites/normalization_check.txt` | `python -m thesis.analysis.check_normalization` |
| 21 | **Differential sites**: 2,973 tested, **1,548 significant** (BH-FDR<0.05) under a batch-aware linear model; 801 up, 747 down; 659 on the shortlist. Welch alone gives 939, so **643 sites (21.6%) change significance** when batch is modelled | `out/diffsites/diffsites.txt` | `python -m thesis.analysis.differential_sites` |
| 22 | **KSEA and ranking are NOT meaningfully disjoint.** Overlap 0/1/3 at the top 10/20/50 cutoffs against 0.56/1.50/3.56 expected by chance; hypergeometric depletion n.s. in all 6 comparisons | `out/ksea_ranking_overlap.txt` | `python -m thesis.analysis.ksea_ranking_overlap` |
| 23 | **Negatives are an alphabetical head slice, not a sample.** Only **55 of 420** kinases are ever a negative; **361 of 415 scored (87.0%) never are**. No significant rank consequence (p 0.060 to 0.840 across 6 tests) | `out/negative_sampling_audit.txt`, `.csv` | `python -m thesis.analysis.negative_sampling_audit` |
| 24 | **Differential sites cannot be protein-corrected.** ProteinExpression is 91.7% empty, 2 rows with >=6 in both groups; **0 of 659 shortlist sites** have usable protein data. Presence/absence sites are filter artifacts but NOT batch artifacts (0/114 confined to one batch, 114/114 had the opposite group covered) | `out/diffsites/diffsites_limitations.txt` | `python -m thesis.analysis.diffsites_limitations` |
| 25 | **Random negatives: the GLOBAL effect is real.** 420/420 kinases used as negatives vs 55. Generic moves 79->74 (r=+0.68, stable over 3 sampler seeds); Liver FC barely moves. Context effect shrinks +0.380 -> +0.283 | `out/random_negatives/randneg_compare.txt`, `randneg_vs_original.csv`, `randneg_context_effect.csv` | `python -m thesis.analysis.random_negatives --mode loo --network <net> --seed <n>` then `python -m thesis.analysis.randneg_compare` |
| 26 | **Fixing the sampler makes site-invariance WORSE.** eta^2 0.9997 -> **1.0000**, median per-kinase rank SD 1.05 -> 0.16, kinases with SD<=1 48.2% -> 93.3%. Popularity baseline global advantage unchanged (-0.758 -> -0.742) | `out/random_negatives/randneg_site_dependence.csv`, `randneg_vs_popularity.csv` | same |

| 27 | **The spectral embedding ships sqrt(degree) as 12-17 of its 32 dimensions.** `SpectralEmbeddingStrategy` requests `dims+1` eigenvectors and drops exactly ONE, correct only for a connected graph. The embedding graphs have **34 / 116 / 122** components of size >=2 (generic / control / liver FC), so the normalized Laplacian's zero eigenvalue has that multiplicity: **13 / 18 / 17** of the 33 computed eigenvalues are numerically zero (~1.7e-08) and the first genuinely non-zero one is 3.53e-02 / 3.65e-03 / 3.68e-03. Null vectors restricted to the giant component are scalar multiples of sqrt(deg): **R^2 = 1.0000** for the shipped null dims, against 0.0011 for generic's informative dims; control's **median over all 32** shipped dims is R^2 = 1.0000. Keeping only eigenvectors with eigenvalue > 1e-6 gives **0/32** degree dims | `out/relational_only/spectral_spectrum.txt`, `spectral_nullspace.txt`, `spectral_fix.txt` | `python -m thesis.analysis.spectral_spectrum`, `.spectral_nullspace`, `.spectral_fix` |
| 28 | **The predicted ranking IS graph degree.** Spearman rho between a kinase's embedding-graph degree and its mean `rank_in_site` in the existing Phase 3 output: **-0.9995** (control), **-0.9973** (liver FC), -0.9387 (generic); Pearson R^2 0.97 / 0.97 / 0.85. PSP substrate count, the popularity baseline's key, reaches only rho = -0.44, R^2 = 0.13. Degree is computed on the embedding graph, which excludes `phosphorylates`, so this is not label leakage. This is a single mechanical cause for findings 9, 13, 14 and 26 and for the Phase 2 hub bias | `out/relational_only/degree_rank.txt`, `degree_rank_{network}.csv` | `python -m thesis.analysis.degree_rank` |
| 29 | **Site embeddings are collapsed under the context networks, and removing the null space does not fix it.** Share of site-embedding magnitude common to ALL 3,228 sites: **0.9160** (liver FC) against **0.1239** (generic, corrected). **95%** of the Hadamard `K*S` feature vector is a per-kinase constant no weight vector can make site-dependent (median over kinases; control 0.88, generic 0.47). Survives both fixes: dropping all zero eigenvectors leaves liver FC at 0.8882, restricting to the giant component at 0.8882. The spectral gap is the difference: 3.68e-03 against generic's 3.53e-02 | `out/relational_only/embedding_geometry.txt`, `spectral_fix.txt`, `spectral_nullspace.txt` | `python -m thesis.analysis.embedding_geometry` |
| 30 | **The degree channel is SPECTRAL-SPECIFIC, and spectral is the least site-aware embedding.** degree-vs-mean-rank rho: spectral -0.9387 / -0.9995 / -0.9973 (generic / control / liver FC) against node2vec **-0.3206 / -0.4362 / -0.3594** and concat -0.3252 / -0.4248 / -0.3824. Median per-kinase rank SD: spectral 4.93 / 1.05 / 1.18 against node2vec **10.20 / 9.08 / 8.20**. Node2vec ranks vary 8-10 positions by site; spectral's vary 1. The prior expectation that node2vec would encode degree comparably was WRONG. **Phase 2's reason for preferring spectral does not survive**: it halved the median rank by encoding a degree prior while being site-blind. **CAVEAT: all rank SDs are single-seed, and findings 11 and 12 put node2vec's seed-to-seed variation above the reported effect sizes, so part of node2vec's 8-10 position spread may be random-walk noise rather than site information. "node2vec is more site-aware" is NOT established.** The spectral side is unaffected (deterministic given the graph), so what stands is that spectral's ranks barely move across sites | `out/relational_only/degree_rank.txt` | `python -m thesis.analysis.degree_rank` |
| 31 | **Kinase degree is 95.5% MSigDB pathway cliques.** Edges touching a candidate kinase, PSP excluded: `protein_same_pathway` 273,279 (**95.5%**), `ppi_high_confidence` 11,115 (3.9%), `has_site` 1,824 (0.6%), `dephosphorylates` **1**. So the ranking is pathway-annotation breadth, which closes the chain to the 90.6% clique finding: cliques -> degree -> null-space dims -> rank | `out/relational_only/degree_composition.csv` | same |
| 32 | **The site collapse is caused by the correlation edges and is reversible; a MAGNITUDE threshold is the right filter, not a sample floor and not FDR.** site_rho / spectral gap / components / sites kept: current 80th pct 0.8882 / 3.68e-03 / 122 / 3,228; FDR<0.05 0.8462 / 3.93e-03 / 122 / 3,228 (removes only 13%, **86.7% of the already-thresholded pool passes**); n_common>=12 0.2244 / 5.53e-03 / 338 / 3,228; **top 25% by |r| 0.3685 / 8.08e-03 / 127 / 3,228**; n_common=18 0.1622 / 1.14e-02 / 667 / 2,529; no corr edges 0.0588 / 3.02e-02 / 998 / **1,860**. A sample floor strands sites detected in few samples; a magnitude threshold prunes uniformly. Confirms CLAUDE.md's base-relations-only warning: 42% of sites stranded | `out/relational_only/corr_edge_variants.txt`, `corr_edge_variants_extra.csv` | `python -m thesis.analysis.corr_edge_variants --k 96` |
| 33 | **The null-space fix works on the embedding and FAILS on the ranking.** Branch `fix/spectral-nullspace`, `drop_null_space=True`, everything else identical (frozen trials, 415 candidates, original sampler, d=32, seed 42). Degree dimensions over giant-component nodes: **17/32 -> 0/32** (control), **16/32 -> 0/32** (liver FC). But degree-vs-mean-rank Spearman moved only **-0.9995 -> -0.9797** (control) and **-0.9973 -> -0.9802** (liver FC). LOO median adjusted rank **unchanged**: generic 79.0 -> 79.0, liver FC 74.0 -> 74.0. Paired Wilcoxon on 6,894 trials favours the ORIGINAL on generic (median diff -2.0, rank-biserial -0.625) and is a tie on liver FC (median diff 0.0, 4,499 of 6,894 trials identical). So removing the dimensions that ARE degree does not remove degree from the prediction | `out/spectral_fixed/spectral_fixed_compare.txt`, `fixed_vs_original.csv`, `fixed_degree_eta.csv` | `python -m thesis.analysis.spectral_fixed_runs --mode loo --network generic` (x4), then `python -m thesis.analysis.spectral_fixed_compare` |
| 34 | **The surviving degree channel is the embedding-vector NORM, and it is intrinsic to the normalized Laplacian.** `L_norm u = lambda u` is equivalent to `L_rw (D^-1/2 u) = lambda (D^-1/2 u)`, so every eigenvector is `u = D^1/2 v` and entry i of ANY of them carries a factor sqrt(deg_i). The null space is only the case where v is constant and u is EXACTLY sqrt(deg). Measured Spearman between a kinase's embedding-vector norm and its degree: **+0.9795 -> +0.9438** (control) and **+0.9444 -> +0.9444** (liver FC) across the fix. Scaling rows to unit length (Ng-Jordan-Weiss) is the candidate remedy and removes the norm by construction; whether it changes held-out ranks is **UNTESTED** | `out/spectral_fixed/whynot_degree_scale.txt`, `whynot_degree_scale.csv` | `python -m thesis.analysis.spectral_fixed_whynot` |
| 35 | **The fix makes site-invariance WORSE, the same direction as the random-negatives fix (finding 26).** eta^2 **0.9997 -> 0.9999** (control) and **0.9997 -> 1.0000** (liver FC); median per-kinase rank SD **1.05 -> 0.52** and **1.18 -> 0.24**; share of kinases with SD <= 1 rank **48.2% -> 76.6%** and **42.2% -> 93.7%**. Two independent, separately-motivated corrections have now both increased site-invariance, which makes a shared downstream cause more likely than a coincidence | `out/spectral_fixed/fixed_degree_eta.csv` | same |
| 36 | **The popularity baseline still beats the fixed models, and the poor-kinase exception survives.** Global: popularity median 25.0 against fixed generic 79.0 (rank-biserial **-0.750**, was -0.758) and fixed liver FC 74.0 (**-0.687**, was -0.686). Unchanged within noise. The exception holds: on **poor** kinases the fixed models win globally, generic **+0.627** (median 303.0 against 233.0) and liver FC **+0.481** (303.0 against 273.5), n = 198. Within-category the baseline wins everywhere, including poor | `out/spectral_fixed/fixed_vs_popularity.csv` | same |
| 37 | **The context effect does not shrink under the fix; it grows slightly.** Generic vs Liver FC on the same 6,894 trials: rank-biserial **+0.380 -> +0.428**, median paired difference **+6 -> +9** ranks, both p < 1e-200. Direction by category is unchanged, Liver FC better on rich (+0.594) and Generic better on average (-0.376) and poor (-0.914). So the context effect is not an artifact of the null-space dimensions specifically. It remains subject to finding 11, seed noise (+0.373) exceeding the effect | `out/spectral_fixed/fixed_context_effect.csv` | same |
Provenance (git hash per results file): `out/stat_tests_provenance.csv`.

## Retracted or reframed claims

1. **Poster/README node2vec numbers** (188 vs 165–174) — single seed; cross-seed spread exceeds the difference.
2. **Hub-monopoly tables** (ULK1 100%, VRK2 94.8%, VRK3 95.5%) — seed 42 only; six runs gave six winners.
3. **Necroptosis cluster as liver-specific** — degree-matched null reproduces 100% of the shift.
4. **"Spectral halves median rank"** — trial-level only; not significant per kinase.
5. **KSEA 17/64** — use **12/64** two-tailed; cite 17/64 once as the original convention.
6. **`batch_words=4` caused the 343.8-min runtime** — false. Correct-settings runs took 357–378 min (3.1–3.3 s/trial); the trial loop dominates.
7. **|r| threshold admits low-n pairs** — false. Mean |r| is flat across n buckets, highest at n=18.
8. **Cluster/all-sites results** — excluded; runs were incomplete (timeouts). Do not cite `allsites_spectral_summary.txt`.
9. **TTK as a surviving exception to the null** — retracted 2026-10-25. Its
   Control to Liver FC shift was **1 rank** (median 6 to 7). A 0% null
   reproduction of a 1-rank shift is noise, not evidence. State it as: TTK's
   shift was too small to evaluate against the null.
10. **Negative sampling: the global effect was tested and is real; the
    direction does not rescue ranking.** RESOLVED 2026-10-26 by findings 25
    and 26. The local effect remains undetected (p 0.060 to 0.840). The
    global effect is measurable: replacing the sampler moves Generic by 5
    rank positions (r = +0.68, stable over 3 sampler seeds) and shrinks the
    context effect by about a quarter. But site-invariance INCREASES
    (eta^2 -> 1.0000) and the popularity baseline's global advantage is
    unchanged, so four results are now known NOT to be sampler artifacts:
    the fixed global ordering, the baseline's strength, the models' edge on
    poor kinases, and the trial-level context effect. Superseded guidance
    retained below for the record:

    **(superseded) Claim neither bias nor absence of bias.** The defect
    is real (finding 23). What has been tested is only the **local** effect,
    whether the 55 ever-negative kinases are pushed down relative to the other
    360, and it was not detected (p 0.060 to 0.840). That comparison **cannot**
    detect the **global** effect: if the classifier only ever separates
    positives from the same 55 early-alphabet kinases, it learns a direction in
    feature space that shifts every kinase's score, and that could itself
    produce the popularity-like prior seen in findings 13 and 14. The global
    effect is **untested** until a rerun with random negatives. Until then
    state: local effect not detected, global effect untested. An earlier note
    here claimed the rank consequence was simply absent; that was too strong.

11. **Misattribution in commit `ac46b90`.** Its message is headed "Correction,
    accepting the advisor's point". That is wrong: the negative-sampling
    correction, both the local-vs-global distinction and the point that a
    rerun changing only the negatives is a valid paired comparison, came from
    the **user**, not from the advisor. History is not rewritten; the
    correction is recorded here and in the message of the commit that follows
    it. Do not cite `ac46b90` as evidence of an advisor request.

12. **`SpectralEmbeddingStrategy`'s docstring is false as implemented.** It claims "High-degree nodes are penalized by construction" and that spectral embeddings do "not share" node2vec's hub bias, and `CLAUDE.md` repeats the second claim. Finding 27 shows 12-17 of the 32 shipped dimensions ARE sqrt(degree), and finding 28 shows the resulting ranking is degree-ordered at rho = -0.9995. The normalization argument is correct for the non-trivial eigenvectors; it does not hold for the null space the code fails to discard. Finding 30 shows both claims are not merely unsupported but **inverted**: spectral is the MORE degree-dependent embedding (-0.9995 against node2vec's -0.4362) and the least site-aware (rank SD 1.05 against 9.08). Do not repeat either claim, and do not cite "spectral fixes the hub bias".

## Quarantined

`out/null/_INVALID_batch_words4/` — node2vec null seed 1 run with `batch_words=4` instead of 10000. Not comparable; see its README. Spectral nulls unaffected.

## Advisor direction (email, 2026-10-15)

Supersedes open decision 1. Every number here was recomputed from the repo on
2026-10-17; commands are read-only.

1. **Journal, not conference.** The grant program manager requires a journal
   paper for the NIH renewal, so there is no deadline pressure. A journal will
   demand the robustness checks already in the findings table, which makes the
   validation suite an asset rather than a liability.

2. **New deliverable: experimental validation targets.** A cancer-center
   collaborator will validate kinases and kinase-substrate pairs, starting with
   substrates that are statistically significant and liver-cancer associated.
   **No per-site differential test exists anywhere in `src/` or
   `src_prediction/`** (grepped `ttest`, `mannwhitney`, `fdr`, `multipletests`,
   `benjamini`: zero hits). This is new work. It is also independent of the
   graph and the model, so it survives every objection the validation round
   raised. The kinase half is partly covered by KSEA (17/64 at adj_p < 0.05
   one-tailed, 12/64 two-tailed), with the caveat that none of those kinases
   reach top-1/5/10 in the LOO rankings.

3. **DEPOD discrepancy confirmed.** `data/raw/PPase_protSubtrates_201903.xls`
   has 938 rows and **115 distinct phosphatases**, matching the advisor's count
   of 114. The 250 previously reported by a prior student is not reproducible
   from this file. Our pipeline is lossier still: 115 to 82 (rows naming a
   substrate but no specific dephosphosite), then 82 to **65 in the built
   graph**, because `_has_exactly_one_site_token` in
   `src/builders/ppase_substrate.py` drops any row listing multiple
   dephosphosites. Blast radius is small: `dephosphorylates` is **329 of
   1,095,542** generic edges (0.03%), so it cannot explain any measured result.
   Worth fixing for correctness, not a reason to redo all preprocessing.

4. **Generic graph composition** (from `outputs/edges.csv.gz`, total 1,095,542).
   This is the more consequential preprocessing question, and it is consistent
   with why the degree-preserving null reproduces our results.

   | Relation | Edges | Share |
   |---|---|---|
   | `protein_same_pathway` | 992,485 | **90.6%** |
   | `ppi_high_confidence` | 36,009 | 3.3% |
   | `site_coevolution` | 24,005 | 2.2% |
   | `site_same_pathway` | 17,922 | 1.6% |
   | `phosphorylates` (target) | 14,596 | 1.3% |
   | `has_site` | 10,196 | 0.9% |
   | `dephosphorylates` | 329 | 0.03% |

5. **Correlation sample counts, the advisor's explicit question.** He asked
   "Don't we calculate correlation across all samples?" The answer is no.
   `src/graph/liver_network.py:323` calls
   `fc_mat.T.corr(method="pearson", min_periods=min_common)`, and pandas uses
   **pairwise complete observations**, so each site pair is correlated using
   only the samples where both sites were quantified. `src/runners/run_liver.py`
   passes only `site_corr_percentile=80.0`, so `min_common` and
   `min_identified` both sit at the default of **6**. Distribution over the
   976,855 stored FC edges: **n=6 is the mode at 151,197 (15.5%)**, n=7-11 is
   445,875 (45.6%), n=12-17 is 301,645 (30.9%), **n=18 is 78,138 (8.0%)**.
   The methodological problem is not the counts. It is that the
   80th-percentile cut is computed over all pairs pooled regardless of n, so
   an r of 0.9 from 6 points competes on the same scale against an r of 0.9
   from 18 points, with no significance control anywhere. Mean |r| is flat
   across n buckets (retraction 7), so the filter does not preferentially
   admit low-n pairs, but it does treat them as equally trustworthy.

6. **Analysis redirection: within-category comparative ranking.** The advisor
   wants prediction performance dropped as the headline metric, replaced by
   comparative kinase ranking with kinases ranked inside their own
   substrate-count category. This is **already implemented** as
   `adjusted_held_out_rank_in_category` (Phase 2, since June), and finding 6
   above is the result: rich ties (33 to 32, not significant), while average
   and poor are **significantly worse** under context. Within-category ranking
   does fix the denominator problem he raised (poor kinases go from a global
   median of 231 to 75), but it removes the only apparent support for the
   hypothesis rather than rescuing it. Separately, comparative ranking is the
   analysis the degree-preserving null reproduces (findings 9, 10 and 13), so
   changing the headline metric does not escape that. Both points need to be
   on the table before the pivot is agreed.

7. **Requested: a step-by-step walkthrough meeting** on graph construction and
   integration, to get aligned. Not yet scheduled, no walkthrough document
   drafted.

## Advisor decisions (2026-10-25)

Received after she read `ADVISOR_BRIEF.md`. These supersede the 2026-10-15
notes above where they conflict.

1. **Journal manuscript, not a conference paper.** The grant program manager
   requires a journal paper for the NIH renewal. No deadline pressure.
2. **Comparative kinase ranking with within-category ranking, not PSP
   prediction performance.** PhosphoSitePlus is not liver-specific, so
   "performance improved with liver data" does not test our hypothesis.
   Kinases are to be ranked inside their own substrate-count category so
   understudied kinases do not compete with rich ones.
3. **Collaborator deliverable.** The cancer center wants statistically
   significant, liver-cancer-associated **substrates** for experimental
   validation, kinases and kinase-substrate pairs both.
4. **Preprocessing is under suspicion.** She found 114 phosphatases in DEPOD
   against the 250 reported earlier and wants the preprocessing re-examined.
   Confirmed: the file holds 115, our graph retains 65. See finding 18 for a
   second, independent preprocessing defect found during the audit.

## Current plan

| # | Task | State | Output |
|---|---|---|---|
| 1 | Record these decisions in STATE.md | **done** 2026-10-25 | this file |
| 2 | Within-category comparative ranking vs the null, read-only | **done** 2026-10-25 | `out/withincat_null.txt` + 3 CSVs, finding 17 |
| 3 | Preprocessing audit written for a biologist | **done** 2026-10-25 | `thesis/GRAPH_WALKTHROUGH.md` |
| 4 | Differential site analysis for the collaborator | **approved and run** 2026-10-25, with per-channel normalization and a batch-aware model added | `thesis/PLAN_differential_sites.md`, `out/diffsites/` |
| 5 | Chapter III updated with the audit; Chapter IV drafted; Chapter II citations and APA fixed | **done** 2026-10-25 | `thesis/chapters/` |
| 6 | Negative-sampling audit; Chapter IV corrections; diffsites limitations; **Chapter V drafted** | **done** 2026-10-26 | `thesis/chapters/chapter_05_*.md` |

Task 2's result, stated plainly: within-category ranking is the right metric and
it does **not** rescue the context effect. Of 415 kinases, 184 (44.3%) have a
Liver FC within-category rank distribution that differs from the
degree-preserving null at BH-FDR < 0.05, but that significance is a sample-size
artifact. Each kinase is tested across 3,228 sites, so the KS test resolves
differences far below one rank position. Among those 184 kinases, 159 have a
median shift of exactly **0** ranks, 23 have 1, one has 2, and none exceeds 2.
Across all 415 kinases the Liver FC shift and the null shift correlate at
**r = 0.9944** (Spearman 0.9600), and the median gap between them is 0.
Separately, Control vs Liver FC is significant for 75.9% of kinases and Control
vs Null for 74.9%, which are the same number within noise.

Task 4 is blocked on approval, and partly on the patient pairing: the design is
consistent with 18 matched pairs but the file does not record which control
channel belongs to which tumour channel.

## Open decisions

1. **MOVED 2026-10-25 to "Advisor decisions" above.** Options A-D in
   `ADVISOR_BRIEF.md` are retired. The direction is settled: journal
   manuscript, within-category comparative ranking as the headline analysis,
   and the collaborator's differential substrates as the new contribution.
   What remains open is narrower and is now decision 8.
2. **PARTIALLY LIFTED, then the hypothesis SUPERSEDED 2026-10-26 by findings
   27-32.** The binding constraint is not the feature layout but the spectral
   embedding's Laplacian null space: the ranking is graph degree at
   rho = -0.9995, and the site embeddings are collapsed (site rho 0.8882), so
   Option D cannot recover site specificity that is absent from the embedding.
   Option D is **not disproved** as an additional contributor, and **no Option D
   run was started**. See `thesis/EMBEDDING_DEFECT.md` and open decision 13.
   `CATCH_UP.md` (c) was rewritten to lead with the degree finding and (e) no
   longer presents Option D as the next experiment. Original text follows.

   The **relational-only baseline** now has a
   source: the user named it **Option D** and asked for it in writing as the
   proposed next experiment, so `CATCH_UP.md` (e) describes it. It drops the K
   and S blocks from `[K, S, |K-S|, K*S]` and keeps only `|K-S|` and `K*S`, so
   no feature is constant across sites for a given kinase. Not yet run; it
   needs a change to `src_prediction/data_prep/pair_features.py`, so it waits
   on approval. The hypothesis it tests, also the user's and recorded in
   `CATCH_UP.md` (c) as a leading hypothesis rather than a fact: the K block
   lets the classifier learn a per-kinase base rate that is added at every
   site, which would explain finding 14 (the popularity baseline's strength)
   and finding 26 (eta^2 rising to 1.0000 once uniform negatives let that base
   rate be estimated cleanly for all 420 kinases). Supporting structural facts:
   `pair_features.py` builds `np.hstack([K, S, abs(K-S), K*S])`, and under the
   random sampler each kinase is a negative 242 to 344 times
   (`out/random_negatives/randneg_loo_generic_spectral_seed2_negative_usage.csv`)
   against 2,458 to 2,459 times for the 55 head-slice kinases before.
   **Still held, no source:** context-conditioned kinase feature, the
   15.2%-vs-60% gate, drop-K and the lbfgs flag.
3. **MOVED 2026-10-25.** The chapter hold is lifted in principle: the
   direction is settled enough to frame Chapters I and V. Chapter VI still
   waits on the differential-site result, since that is now the contribution
   it would report. The walkthrough meeting has not happened yet, but
   `thesis/GRAPH_WALKTHROUGH.md` is written and ready for it.
4. **G5–G14 doc fixes** — RESOLVED 2026-10-17. Committed as `c5f6111` and
   merged into `main` by fast-forward, no merge commit. `main` is now at
   `b72234c`.
5. **Provenance of the liver dataset. Narrowed 2026-10-25, still open.**
   The sheets carry their original captions, which are the strongest lead yet:
   - `ProteinExpression`: "Table S4 Detailed information of unique proteins
     identified in HCC tumour tissues"
   - `Phosphorylation`: "Table S11 Identified and quantified phosphopeptides,
     phosphosites and corresponding phosphoproteins in HCC tumour and adjacent
     non-tumour tissues"
   **Gao et al. (2019) is ruled out**: that cohort is 159 patients, ours is 18
   and 18. Two extended web searches on the captions did not identify the
   paper. Next step is to ask the advisor or the file's supplier, not to search
   again. Until then, describe the data by its contents and sample counts and
   cite nothing. British spelling in the captions suggests a non-US journal.
   **Implemented 2026-10-25**: Chapters II and III now describe the dataset
   without attribution, and `references.md` records the evidence gathered.
6. **RESOLVED 2026-10-25. Two unsupported claims in the chapter drafts.**
   Both were addressed:
   - The 85th-percentile sentence was **deleted** from Chapter III. A repository
     search found no output file for that comparison, only the assertion in
     `CLAUDE.md`, which is itself the unsupported documentation.
   - The KSEA-disjointness sentence was **computed and corrected** (finding 22).
     The sets overlap in 0, 1 and 3 kinases at the top 10, 20 and 50 cutoffs,
     against 0.56, 1.50 and 3.56 expected by chance, with no significant
     depletion in any of the six comparisons. Chapter II now states that the two
     methods nominate few kinases in common but no more rarely than chance, and
     that the apparent disagreement follows from both nominating small sets from
     a large pool.

   Original wording retained for the record:
   - Chapter III §3.6.2 — "An 85th-percentile threshold was tested on the
     fold-change variant and performed worse." Referenced in project docs but
     no output file for that comparison has been found.
   - Chapter II §2.7 — "KSEA and ranking nominate largely disjoint kinase
     sets." Plausible from `outputs_prediction/ksea/ksea_comparison.txt` and the
     LOO top-1 tables, but the disjointness was never computed directly. Either
     compute it after the advisor meeting or soften to what the files state.
7. **Equation rendering.** Chapter II contains LaTeX (`$$L = I - D^{-1/2} A
   D^{-1/2}$$`, the KSEA score). python-docx has no LaTeX support, so plan the
   conversion for Checkpoint 4: either OMML via a converter, or rendered images
   with text fallback. Decide before the build, since it affects numbering and
   captioning.
8. **What the headline analysis now is.** Finding 17 means within-category
   comparative ranking cannot carry the manuscript on its own: it is
   indistinguishable from a degree-preserving null at any effect size that
   matters. Either the graph changes (see decision 9) or the differential-site
   contribution becomes the primary result and the ranking work is reported as
   a negative methodological finding. This is the main thing to settle with the
   advisor.
9. **Preprocessing fixes: which, if any, to make.** Three defects are
   documented in `GRAPH_WALKTHROUGH.md` §6. By blast radius:
   - `protein_same_pathway` is **90.6%** of generic edges, all MSigDB cliques.
     This is the only one that plausibly changes a result.
   - The truncated identifier map costs **1.5%** of PPI candidate rows.
   - The DEPOD phosphatase loss (115 to 65) is **0.03%** of edges.
   Recommend auditing in that order rather than rebuilding everything.
10. **Patient pairing for the liver data.** Needed for a paired differential
    test. See `PLAN_differential_sites.md` §1 for the three ways to resolve it,
    including a read-only empirical fallback.
11. **RESOLVED 2026-10-26. The rerun was approved and run** (findings 25,
    26; Chapter V 5.10). Six runs total: LOO Generic and Liver FC at sampler
    seeds 1, 2 and 3, plus predict-all for Generic, Control and Liver FC at
    seed 1. 10.4 h of compute. What remains open is whether to adopt the
    random sampler as the project default, which would touch project code
    and invalidate nothing but would require rerunning node2vec and concat
    for consistency. Original text:

    **(superseded) Whether to fix the negative sampler and rerun.** Finding 23 documents
    that negatives come from 55 of 420 kinases. A corrected sampler would be a
    one-line change, but rerunning would break comparability with every
    existing run against the frozen trial set, and the measured rank
    consequence of the defect is nil. **Advisor decision. No rerun has been
    attempted.** Note the fix would touch project code, which is currently
    out of bounds.
12. **Protein-abundance correction for the differential sites.** Not possible
    with the data in hand (finding 24). Either obtain a matched protein
    quantification matrix for the same 36 channels, or state in the manuscript
    that the reported changes are phosphosite-level and uncorrected.

13. **The embedding defect: (a) is DONE and did not work.** See
    `thesis/EMBEDDING_DEFECT.md` and findings 27-37. Status of the four
    candidates:
    - **(a) DONE, negative result.** The null-space fix is on branch
      `fix/spectral-nullspace` (opt-in `drop_null_space`, default unchanged).
      It removes the degree dimensions (17/32 -> 0/32) and does NOT remove the
      degree-ordered ranking (rho -0.9995 -> -0.9797), does not improve the
      median rank (79.0 and 74.0, unchanged), and makes site-invariance worse
      (eta^2 -> 1.0000). Findings 33 and 35. **Not merged.** The branch is
      worth keeping because the bug is real and the default is safe, but it is
      not the lever.
    - **(a2) NEW candidate, untested: row-normalize the embedding.** Finding 34
      shows the surviving channel is the embedding-vector norm, which carries
      sqrt(degree) for every eigenvector of the normalized Laplacian, not just
      the null space. Unit-length rows remove it by construction. This is the
      natural successor to (a) and is the same size of change. Whether it
      changes held-out ranks is untested, about 1.6 h per LOO run at the
      measured rate.
    - **(b) STILL OPEN AND NOW MORE IMPORTANT.** Correlation-edge threshold
      from the pooled 80th percentile to ~95th. Finding 29's site collapse was
      never attributable to the null space, and (a) failing does not touch it.
      This is a change to the graph build in `src/`.
    - **(c) Redo the spectral-vs-node2vec comparison.** Now motivated
      differently: not because the defect explains spectral's advantage, but
      because the defect is NOT what explains it. Subject to the seed caveat on
      finding 30.
    - **(d) Report degree as a baseline.** Unaffected and still worth doing:
      degree reproduces the ranking at rho = -0.98 even after the fix.
    The ordering claim in the previous version of this entry, that anything
    rerun before (a) and (b) just re-measures the degree channel, was **wrong
    about (a)**: the fixed runs still measure the degree channel, because the
    channel is not in the dimensions that were removed.

## Outstanding work

- Chapter II: citations verified, APA 7 applied, disjointness claim corrected
- Chapter III: audit findings integrated, unsupported sentence removed
- Chapter IV: drafted 2026-10-25, corrected 2026-10-26 (negative sampling in
  4.5 and 4.6.2; seed variance scope in 4.9.2; node2vec walk parameters)
- Chapter V: drafted 2026-10-26, section 5.10 (random negatives) added the
  same day, summary renumbered to 5.11 with 7 new rows. Every number carries its `out/` path in an
  HTML comment; effect sizes accompany every p-value; figures and tables
  proposed per section
- **Chapters I and VI: on hold until 2026-10-31**, pending open decision 8
- Chapters I, IV, V, VI, VII after the direction settles
- `.docx` build + PDF export + formatting check (Checkpoint 4)
- Graph-construction walkthrough: **written**, `thesis/GRAPH_WALKTHROUGH.md`
- Per-site differential analysis: **run 2026-10-26** (findings 20, 21, 24);
  plan in `thesis/PLAN_differential_sites.md`, output in `out/diffsites/`
- Catch-up briefing: `thesis/CATCH_UP.md`, rewritten 2026-10-26 to add the
  K-block mechanism as a leading hypothesis (section c) and Option D as the
  proposed next experiment (section e). ~1,290 words
- Advisor package: `thesis/advisor_package/` holds the 1-page differential
  summary, `differential_shortlist.csv` (659 sites) and a copy of
  `GRAPH_WALKTHROUGH.md`
- **Experiments halted 2026-10-26 at the user's instruction** until after the
  advisor meeting. No new analyses or runs, Option D included
- Decide the headline analysis (open decision 8)
- Preprocessing fixes, if any (open decision 9)

## Not to be rerun

`run_overnight.ps1 -Steps c` (node2vec nulls) remains available but is not
required for any current finding. No experiments should run before the advisor
meeting.
