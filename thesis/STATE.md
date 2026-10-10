# Thesis state — recovery point

Last updated: 2026-10-10. Analysis phase **complete**; no experiments pending.
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

## Open decisions

1. **Advisor meeting** — Option A (evaluation thesis, 0 h), B (all variants, ~14–15 h), C (variant d only, ~2.1 h) or D (relational-only features, uncosted). See `ADVISOR_BRIEF.md`.
2. **Context-conditioned kinase feature / 15.2%-vs-60% gate / relational-only baseline / drop-K / lbfgs flag** — not in the repo. Do not write about these until a source is confirmed.
3. **Chapters I, V, VI** — on hold until after the advisor meeting.
4. **G5–G14 doc fixes** — staged on a branch, not committed; diff to be reviewed.

## Outstanding work

- Chapters II and III drafts (in progress; independent of the above)
- G5–G14 doc fixes branch (in progress)
- Chapters I, IV, V, VI, VII after the advisor decision
- `.docx` build + PDF export + formatting check (Checkpoint 4)

## Not to be rerun

`run_overnight.ps1 -Steps c` (node2vec nulls) remains available but is not
required for any current finding. No experiments should run before the advisor
meeting.
