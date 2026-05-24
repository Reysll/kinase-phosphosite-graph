---
name: check-results
description: Read and summarize current experiment results from outputs_prediction/. Surfaces top-line metrics (mean/median rank, top-k accuracy, per-trial deltas) for any experiment — no manual file navigation needed. Invoke as /check-results [experiment] where experiment is: generic, liver, generic-multi, liver-multi, compare, or "all".
---

The user wants a summary of current experiment results.
$ARGUMENTS specifies which experiment(s) to report on. Default to "all" if empty.

## Result file map

| Experiment | Summary file | Metrics file |
|---|---|---|
| `generic` | `outputs_prediction/generic_trained_model_debug/summary.txt` | `outputs_prediction/generic_trained_model_debug/metrics.csv.gz` |
| `liver` | `outputs_prediction/liver_trained_model_fc_corr_debug/summary.txt` | `outputs_prediction/liver_trained_model_fc_corr_debug/metrics.csv.gz` |
| `generic-multi` | `outputs_prediction/generic_multi_kinase/summary.txt` | `outputs_prediction/generic_multi_kinase/metrics.csv.gz` |
| `liver-multi` | `outputs_prediction/liver_multi_kinase/summary.txt` | `outputs_prediction/liver_multi_kinase/metrics.csv.gz` |
| `compare` | `outputs_prediction/comparison_summary_fc_corr.txt` | `outputs_prediction/comparison_improved_in_liver_fc_corr.csv.gz` |

## Instructions

1. Parse $ARGUMENTS to decide which experiment(s) to report. If empty or "all", report every available one.
2. For each requested experiment:
   a. Check if the summary file exists. If not, tell the user the experiment hasn't been run yet.
   b. Read the summary .txt file and extract: number of LOO trials, mean held-out rank, median held-out rank, top-1 / top-5 / top-10 / top-20 accuracy, best-true rank if available.
   c. If a metrics.csv.gz exists, read it to surface per-trial distribution if the summary doesn't already include it (min, max, std of held-out rank).
3. For `compare` or `all`: also read `comparison_summary_fc_corr.txt` and report:
   - How many trials changed top-1 prediction (generic → liver)
   - How many trials had improved (lower) held-out rank in liver
   - The net delta in mean and median rank
4. Present results in a clean markdown table. Group by experiment. Highlight the single most notable finding in bold.
5. At the end, note which files were read and their modification timestamps so the user knows how fresh the numbers are.

## Context for interpretation

- Candidate kinase pool: 420 kinases — random top-1 baseline ≈ 0.24%, top-20 ≈ 4.8%.
- 50-trial single-site LOO: liver (mean 118.2, median 104.0) vs generic (mean 157.4, median 118.5).
- Multi-kinase LOO (6,909 trials): generic tends to win on hub sites (AKT1-S473, GSK3B-S9); liver wins on sparse sites.
- TTK has 68 substrates in the graph and is a liver cancer kinase of particular interest.
