# Plan only — correlation-edge improvement experiment (spectral, null-controlled)

Status: **not run.** No code or data outside `thesis/` would be modified.

## Motivation

Two results from the current evidence set drive this plan.

1. The degree-preserving null reproduces the Liver FC advantage over Generic
   almost exactly under spectral (ranks identical in ~96% of trials), and
   reproduces **100%** of the necroptosis cluster's Control→Liver FC rank shift.
   The current correlation edges therefore act mainly through *how many* edges
   each site gains, not *which* sites are connected.
2. The current edge set is built from very small samples. Of the 976,855 FC
   correlation edges, only **78,138 (8.0%)** rest on the full n = 18; **898,717
   (92.0%)** use fewer, with the most common value at the `min_periods` floor of
   **n = 6** (151,197 edges). A Pearson r estimated on 6 points has very high
   variance, so a large fraction of these edges may be noise.

   **Correction to an earlier claim.** I previously suggested the
   |r|-percentile threshold would preferentially admit low-n pairs, because
   small samples reach extreme |r| more easily. Measured on the retained edges,
   that is **not** supported: mean |r| is essentially flat across sample-size
   buckets and is in fact highest at n = 18.

   | n_common | edges | mean abs(r) |
   |---|---|---|
   | 6-7 | 256,805 | 0.8555 |
   | 8-11 | 340,267 | 0.8370 |
   | 12-14 | 174,904 | 0.8368 |
   | 15-17 | 126,741 | 0.8576 |
   | 18 | 78,138 | 0.8750 |

   The remaining concern is therefore estimator *variance* at small n, not
   systematic inflation of |r|. Note this table is conditional on edges that
   already passed the 80th-percentile cut; the unconditional distribution was
   not measured.

3. Kinase identity alone explains **99.97%** of the variance in `rank_in_site`
   under spectral (eta^2 = 0.9997 for Control and Liver FC; median per-kinase
   rank SD across all 3,228 sites is 1.05-1.18 positions). The scorer applies a
   near-fixed global ordering of kinases and barely conditions on the site.
   Any edge-set variant that leaves this property intact cannot produce
   site-specific predictions, regardless of how the edges are filtered.

Every variant below is therefore tested **against its own degree-preserving
null**, because a variant that only changes edge count will be indistinguishable
from its null.

## Variants

All counts computed from the stored `weight` (Pearson r) and `n_common` columns
of `outputs/Liver_network_edges.csv.gz`; no recomputation was needed.

### (a) FDR-corrected significance of each correlation

Per-pair two-sided t-test on r with that pair's own `n_common` (not a blanket
n = 18), Benjamini-Hochberg across all stored FC pairs.

| Threshold | Edges kept | % of current |
|---|---|---|
| raw p < 0.05 | 857,764 | 87.8% |
| FDR < 0.05 | **846,586** | 86.7% |
| FDR < 0.01 | 673,519 | 68.9% |

Caveat: these are computed over the already-thresholded top-20% set, so they are
conditional on the existing filter. A clean version requires recomputing
correlations for all pairs from the Excel, which is the expensive path below.

### (b) Stricter |r| percentile

Current threshold is the 80th percentile (top 20%). Subsetting the stored edges
by |r| gives the stricter cuts directly.

| Variant | |r| threshold | Edges kept | % of current |
|---|---|---|---|
| top 5% overall | ≥ 0.9336 | **244,214** | 25.0% |
| top 1% overall | ≥ 0.9797 | **48,843** | 5.0% |

### (d) Minimum overlapping-sample count

Keep only edges whose correlation rests on enough samples to be estimable,
alone and intersected with the top 5% / top 1% |r| cuts from (b).

| Variant | Edges kept | % of current |
|---|---|---|
| n_common >= 12 | **379,783** | 38.9% |
| n_common >= 15 | **204,879** | 21.0% |
| n_common == 18 (complete only) | 78,138 | 8.0% |
| n >= 12 AND top 5% |r| | **77,356** | 7.9% |
| n >= 15 AND top 5% |r| | 40,925 | 4.2% |
| n >= 12 AND top 1% |r| | 9,492 | 1.0% |
| n >= 15 AND top 1% |r| | 4,439 | 0.5% |

Recommended pair to test first: `n >= 12` alone (keeps 38.9%, the mildest
intervention) and `n >= 12 AND top 5%` (keeps 7.9%, the aggressive one). Each
against its own degree-preserving null, same success criterion.

Runtime: 4 additional spectral LOO runs, ~62 min each, ~4.2 h.

### (c) Differential correlation, tumor vs control

Fisher z-test on (r_cancer − r_control) per site pair, FDR-corrected. This is
the only variant that targets *disease-specific rewiring* rather than
correlation strength.

**Edge count cannot be derived from stored data.** Stored cancer and control
edges are each already top-20% filtered, so their intersection is biased and the
pairs that are strong in only one condition — exactly the differential ones —
are partly missing. This variant requires recomputing both full correlation
matrices from the proteomics Excel.

## Code touched

| Item | Reuse as-is | New in `thesis/` |
|---|---|---|
| Edge filtering (a), (b) | reads `outputs/Liver_network_edges.csv.gz` read-only | `thesis/analysis/edge_variants.py` |
| Differential corr (c) | `src/graph/liver_network.py` correlation helpers, imported read-only | `thesis/analysis/diff_corr_edges.py` |
| Null generation | `thesis.analysis.null_network.degree_preserving_swap` | reused |
| LOO | `src_prediction.engine.leave_one_out.run_leave_one_out` | reused unmodified |
| Comparison | `thesis.analysis.null_compare`, `stat_tests` | extended |

Nothing in `src/`, `src_prediction/`, `data/` or `outputs/` is written.

## Runtime

Spectral LOO measured at **~62 min** per 6,909-trial run; null graph build ~95 s.

| Step | Runs | Estimate |
|---|---|---|
| (a) FDR 0.05 + its null | 2 | ~2.1 h |
| (b) top 5% + null, top 1% + null | 4 | ~4.2 h |
| (c) differential + null (incl. ~1–2 h recomputation) | 2 + prep | ~3.5–4.5 h |
| **Total** | 8 LOO runs | **~10–11 h** |

Sparser graphs should run somewhat faster than the measured 62 min, so this is
an upper bound. node2vec is deliberately excluded: at ~6 h per run it would add
~48 h, and the seed-variance result shows single-seed node2vec cannot support
these comparisons anyway.

## New files written (all under `thesis/`)

```
thesis/analysis/edge_variants.py
thesis/analysis/diff_corr_edges.py
thesis/analysis/out/variants/{variant}_edges_meta.json
thesis/analysis/out/variants/{variant}_results.csv.gz
thesis/analysis/out/variants/{variant}_null_results.csv.gz
thesis/analysis/out/variants_compare.{txt,csv}
```

## Success criterion, stated in advance

A variant counts as an improvement only if it clears **both** bars. This applies
to every variant (a)-(d) and to Option D in the advisor brief.

**Bar 1 — its own degree-preserving null.** The variant must beat a null built
from its own edge set with degrees preserved, by more than the seed-variance
ceiling for that embedding. For spectral that ceiling is near zero
(deterministic), so the test is whether the variant-vs-null contrast rises above
the ~96%-tie regime seen with the current edges. Beating Generic is not
evidence, since the current null already does that with the same effect size
(+0.380).

**Bar 2 — the popularity baseline** (`thesis/analysis/out/popularity_baseline.txt`,
parity-corrected: substrate counts from the model's own 6,909 training edges,
same 415 candidates). The variant must beat it:

| Comparison | Popularity | Must beat |
|---|---|---|
| All trials, global rank | median 25.0, Hit@10 33.4%, Hit@20 45.1% | best current model: median 74.0 |
| Poor stratum, global | median 303.0 | models already win here (r +0.48 to +0.81) |
| Poor stratum, within-category | median 64.0 | models currently lose badly (r -0.84 to -0.95) |
| Average stratum, global | median 162.0 | models currently lose (r -0.24 to -0.32) |
| Average stratum, within-category | median 53.0 | models currently lose (r -0.85 to -0.95) |

The poor/average strata are the decisive ones: the all-trials number is
dominated by rich kinases (5,872 of 6,894 trials, 85.2%), so a variant can move
it without improving anything about annotation bias. A variant that beats
popularity **only** on all-trials and not in the sparse strata should be
reported as not having addressed the problem.

## Flags

- No project code or data would be modified.
- Variant (c) is the only one that cannot be costed exactly without the
  recomputation step.
- The n_common finding above is independent of these variants and arguably
  warrants reporting on its own, since it affects every existing correlation
  result in the thesis.
