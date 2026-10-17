# Thesis state — recovery point

Last updated: 2026-10-17. Analysis phase **complete**; no experiments pending.
Advisor response received 2026-10-15; see "Advisor direction" below.
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
| 10 | **Null reproduces 100% of necroptosis shift**; TTK 0% | `out/null_casestudy.txt`, `null_casestudy_ks.csv` | `python -m thesis.analysis.null_casestudy` |
| 11 | Seed noise (+0.373) exceeds context effect (+0.249); seed 13 reverses direction | `out/seedvar_compare.txt` | `python -m thesis.analysis.seedvar_compare` |
| 12 | Six runs, six different top-1 hubs; Generic median spread 40 ranks | `out/seed_averaged.txt`, `seed_averaged_top1.csv` | `python -m thesis.analysis.seed_averaged` |
| 13 | **η² = 0.9997** (spectral); median per-kinase rank SD ≈ 1 across 3,228 sites | `out/site_dependence.txt`, `site_dependence.csv` | `python -m thesis.analysis.site_dependence` |
| 14 | **Popularity baseline beats all models** (median 25.0 vs 74.0; Hit@20 45.1% vs 22.1%) | `out/popularity_baseline.txt` | `python -m thesis.analysis.popularity_baseline` |
| 15 | **Models beat popularity for poor kinases, global rank** (r +0.48 to +0.81) | `out/popularity_parity_stratified.csv` | same |
| 16 | Liver-measured subset (n=510): spectral context survives; null contrast n.s. | `out/subset_analysis.txt` | `python -m thesis.analysis.subset_analysis` |

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

## Open decisions

1. **Direction after the advisor response.** Options A-D in `ADVISOR_BRIEF.md`
   are superseded. The advisor's direction maps to none of them cleanly: it is
   closest to A (no new runs) with within-category comparative ranking as the
   headline, but that headline needs the null and eta^2 results confronted
   first. Pending: whether to scope the differential substrate analysis
   (item 2 above) as the primary new contribution.
2. **Context-conditioned kinase feature / 15.2%-vs-60% gate / relational-only baseline / drop-K / lbfgs flag** — not in the repo. Do not write about these until a source is confirmed.
3. **Chapters I, V, VI** — still on hold. The advisor meeting happened by
   email, but the direction is not settled and the walkthrough meeting
   (item 7 above) has not occurred, so the framing these chapters need is
   still open.
4. **G5–G14 doc fixes** — RESOLVED 2026-10-17. Committed as `c5f6111` and
   merged into `main` by fast-forward, no merge commit. `main` is now at
   `b72234c`.
5. **Provenance of the liver dataset.** Confirm the source study behind
   `data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx` before citing Gao
   et al. (2019). The file has **18 tumor + 18 control** samples, which may not
   match that study's cohort size. If it cannot be confirmed, describe the data
   by its contents and sample counts and drop the attribution.
6. **Two unsupported claims in the chapter drafts.** Locate the file that
   supports each, or delete the sentence:
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

## Outstanding work

- Chapters II and III drafts (in progress; independent of the above)
- Chapters I, IV, V, VI, VII after the direction settles
- `.docx` build + PDF export + formatting check (Checkpoint 4)
- Reply to the advisor's 2026-10-15 email (not drafted)
- Graph-construction walkthrough for the requested meeting (not drafted)
- Per-site differential analysis for the collaborator's validation targets
  (not started; no such code exists yet)
- Decide whether to fix the phosphatase builder's multi-site row loss (65 of
  115 phosphatases retained; 0.03% of edges, so low priority)

## Not to be rerun

`run_overnight.ps1 -Steps c` (node2vec nulls) remains available but is not
required for any current finding. No experiments should run before the advisor
meeting.
