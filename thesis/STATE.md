# Thesis state — recovery point

Last updated: 2026-10-25 (second pass). Advisor decisions recorded below;
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
2. **Context-conditioned kinase feature / 15.2%-vs-60% gate / relational-only baseline / drop-K / lbfgs flag** — not in the repo. Do not write about these until a source is confirmed.
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

## Outstanding work

- Chapter II: citations verified, APA 7 applied, disjointness claim corrected
- Chapter III: audit findings integrated, unsupported sentence removed
- Chapter IV: drafted 2026-10-25, methods only, no results
- Chapters I, IV, V, VI, VII after the direction settles
- `.docx` build + PDF export + formatting check (Checkpoint 4)
- Graph-construction walkthrough: **written**, `thesis/GRAPH_WALKTHROUGH.md`
- Per-site differential analysis: **planned**, `thesis/PLAN_differential_sites.md`,
  awaiting approval before any run
- Decide the headline analysis (open decision 8)
- Preprocessing fixes, if any (open decision 9)

## Not to be rerun

`run_overnight.ps1 -Steps c` (node2vec nulls) remains available but is not
required for any current finding. No experiments should run before the advisor
meeting.
