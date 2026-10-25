# Plan: differential phosphosite analysis for experimental validation

**Status: plan only. Nothing has been run. Awaiting approval.**

Purpose: give the cancer-center collaborator a ranked, statistically defensible
list of phosphosites that differ between HCC tumour and adjacent non-tumour
tissue, with the kinases implicated for each, so they can choose validation
targets.

Why this deliverable is worth doing first: it depends only on the liver
phosphoproteomics file. It does not use the graph, the embeddings or the
ranking model, so none of the problems the validation round found (the
degree-preserving null reproducing our results, seed instability, kinase
identity explaining 99.97% of rank variance) apply to it.

All counts below were measured read-only on 2026-10-25 from
`data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx`.

---

## 1. Are the 18 and 18 matched pairs?

**Checked. Almost certainly yes by design, but the pairing is not recorded in
the file, so a paired test cannot currently be executed.**

The `Phosphorylation` sheet's sample columns, in file order:

```
Control:  F1  F2  F3   F7  F8  F9   F13 F14 F15   F19 F20 F21   F31 F32 F33   F43 F44 F45
Tumour:   F4  F5  F6   F10 F11 F12  F16 F17 F18   F22 F23 F24   F34 F35 F36   F46 F47 F48
```

Six blocks of three controls then three tumours. That is the layout of a
multiplexed experiment with six batches of six channels, three tumour and three
adjacent-normal each, which implies 18 patients each contributing both tissue
types. The sheet caption agrees: "HCC tumour and adjacent non-tumour tissues".

What is missing is the within-batch correspondence. Nothing says whether `F1`
pairs with `F4`, `F5` or `F6`. The `ProteinExpression` sheet uses a completely
different numbering (`F26` to `F61`), so it cannot be used to cross-walk either.

**Three ways to resolve it, in order of preference:**

1. Obtain the sample-to-patient table from the source publication. This is
   blocked on the same question as the citation: we still do not know which
   paper Tables S4 and S11 come from. See `thesis/GRAPH_WALKTHROUGH.md` §5.1.
2. Ask whoever supplied the file.
3. **Empirical fallback.** Correlate every control channel against every tumour
   channel across all sites. If the design is paired, each control should
   correlate most strongly with one specific tumour channel inside its own
   batch, giving a recoverable one-to-one assignment. This is cheap and
   read-only. It would be reported as *inferred* pairing and used only for a
   sensitivity analysis, never as the primary test.

**Recommendation: run the primary analysis unpaired, and treat the paired
version as a sensitivity analysis if and when the pairing is established.** An
unpaired test on a genuinely paired design is conservative, it loses power but
does not inflate the false positive rate. Guessing the pairing wrongly would do
the opposite.

---

## 2. What is actually testable

| Quantity | Count |
|---|---|
| Rows in the `Phosphorylation` sheet | 6,571 |
| Rows parsing to exactly one phosphosite | 5,249 |
| Rows naming more than one site (pipeline drops these) | 1,319 |
| Rows with no parseable site | 3 |
| Overall missingness across the 36 channels | **34.7%** |
| Median observed values per row, control | 13 of 18 |
| Median observed values per row, tumour | 13 of 18 |

Single-site rows passing a minimum-observed filter in **both** groups:

| Minimum per group | Rows |
|---|---|
| 3 | 4,529 |
| 4 | 4,283 |
| **6 (recommended)** | **3,932** |
| 9 | 3,231 |
| 12 (sensitivity tier) | 2,548 |
| 15 | 1,823 |
| 18 (complete only) | 894 |

At the recommended threshold: **2,965 unique phosphosites on 1,682 genes**,
from 3,932 rows. Rows collapse to fewer sites because the same site can appear
on several peptides; see §4 for how duplicates are handled.

---

## 3. Proposed test

### Primary

**Welch's unpaired t-test on log2 abundance**, per site, with Benjamini-Hochberg
FDR across all tested sites.

Reasoning:

- **Log2 first.** Abundances are right-skewed and multiplicative. Everything
  downstream is done on log2 values.
- **Welch rather than Student.** Tumour and adjacent-normal phosphosite
  variances differ, and with 18 per group the unequal-variance correction costs
  almost nothing.
- **Unpaired rather than paired**, for now, for the reason in §1.
- **Parametric rather than rank-based as the primary.** With n up to 18 per
  group the t-test has materially more power than Mann-Whitney, and log2
  abundances are close enough to normal for it. But see the sensitivity test.

### Sensitivity tests, reported alongside, not instead

1. **Mann-Whitney U** on the same sites. Any site significant under the t-test
   but not the rank test should be flagged, because it is probably driven by one
   extreme channel.
2. **Paired t-test and Wilcoxon signed-rank**, only if the pairing is
   established or inferred per §1.3.
3. **The stricter 12-per-group tier** (2,548 sites). A site that survives both
   tiers is a stronger validation candidate.

### On a moderated t-test

A moderated test (empirical Bayes variance shrinkage, as in limma) would be the
best choice statistically, and it is standard in proteomics. I am **not**
proposing it as the primary, for one practical reason: there is no vetted
implementation in this project's dependency set. The options are to add `rpy2`
plus R's limma, which is a heavy new dependency for one analysis, or to write
the empirical-Bayes shrinkage ourselves, which is doable but is new statistical
code that would itself need validating before the collaborator relies on it.

**Recommendation: ship Welch plus BH now, and add a moderated test as a second
pass if the advisor wants the extra power.** I will flag in the output how many
sites sit near the significance boundary, since that is the population a
moderated test would move.

---

## 4. Missing-value handling

**No imputation for the primary test.** Each site is tested on the channels
where it was actually observed.

Phosphoproteomics missingness is mostly not random. A site is usually absent
because its abundance fell below the detection limit, which is informative, so
mean or median imputation would pull low-abundance tumour values up toward the
control mean and bias the result toward the null. Left-censored imputation
exists but it inserts invented numbers into a list a wet lab will spend money
acting on. For this deliverable the honest choice is to test only what was
measured, and to report the observed counts per site in the output so the
collaborator can judge each one.

**Duplicate sites.** Where several peptide rows map to the same site, collapse
by taking the mean of the log2 values per channel, matching the row-mean
collapse `src/graph/liver_network.py` already uses. The collapse happens before
testing, so each site gets exactly one test and the FDR denominator is the
number of sites, not the number of peptides.

**Presence and absence, reported separately.** Sites observed in one group and
never in the other cannot be t-tested, and they are often the most interesting
for validation. These get their own list with a Fisher exact test on the
detection counts, clearly labelled as a detection-rate result rather than an
abundance result.

---

## 5. Effect size

A p-value is not enough to pick validation targets. Each site reports:

| Measure | Why |
|---|---|
| **log2 fold change**, mean(tumour log2) minus mean(control log2) | the interpretable magnitude; the collaborator will filter on this |
| **Hedges' g** | standardized effect, comparable across sites with different observed counts, with the small-sample correction |
| **Cliff's delta** | non-parametric companion, robust to the one-extreme-channel case |
| **n observed per group** | so a site tested on 6 and 6 is never mistaken for one tested on 18 and 18 |

Proposed shortlist rule for the collaborator, to be confirmed with them:
BH-FDR < 0.05 **and** |log2FC| >= 1 **and** at least 9 observed values in both
groups. I will report how many sites survive each clause separately so the
thresholds can be renegotiated without rerunning.

---

## 6. Mapping significant sites to kinases

Two independent routes, deliberately kept separate because they answer different
questions and we already know they disagree.

### Route A: direct PhosphoSitePlus annotation

For each significant site, look up `phosphorylates` edges in
`outputs/edges.csv.gz` to get the kinases annotated to act on exactly that site.

This is the route the collaborator asked for, since they want to validate the
kinase and the kinase-substrate pair together. Expected limitation: PSP
annotates only a minority of measured liver sites, so many significant sites
will have no annotated kinase. The output will state the coverage explicitly
rather than quietly dropping unannotated sites.

### Route B: KSEA on the differential result

Feed the per-site log2 fold changes into the KSEA score already implemented in
`src_prediction/analysis/ksea.py`:

```
score = (mean(substrate log2FC) - mean(all log2FC)) * sqrt(m) / delta
```

This gives per-kinase activity inference across the whole measured set rather
than site by site. Two requirements:

- Use the **two-tailed** p-value, `2 * norm.sf(abs(score))`. The published
  17-of-64 figure used a one-tailed `norm.sf(abs(score))`; the correct
  two-tailed version gives 12 of 64. This is retraction 5 in `STATE.md`.
- The existing module will be **imported and called, never rerun through
  `run_ksea.py`**, because that runner writes into `outputs_prediction/ksea/`.
  All output goes under `thesis/`.

### What will be reported about the disagreement

The two routes nominate different kinases, and `STATE.md` open item 6 notes that
the degree of disjointness was asserted in a chapter draft but never actually
computed. This analysis is the natural place to compute it properly: the overlap
between KSEA-significant kinases and kinases annotated to significant sites gets
reported as a number, which either supports or retires that claim.

---

## 7. Files this would create

All under `thesis/`. Nothing in `src/`, `src_prediction/`, `data/`, `outputs/`
or `outputs_prediction/` is read other than read-only, and none is written.

| Path | Contents |
|---|---|
| `thesis/analysis/differential_sites.py` | the analysis, read-only, runnable as `python -m thesis.analysis.differential_sites` |
| `thesis/analysis/out/diffsites/diffsites.txt` | human-readable report, counts and method record |
| `thesis/analysis/out/diffsites/diffsites_tests.csv` | every tested site: both groups' n, means, log2FC, Hedges' g, Cliff's delta, t and Mann-Whitney p, BH-FDR |
| `thesis/analysis/out/diffsites/diffsites_significant.csv` | the shortlist, sorted for the collaborator |
| `thesis/analysis/out/diffsites/diffsites_presence_absence.csv` | one-group-only sites with Fisher exact results |
| `thesis/analysis/out/diffsites/diffsites_kinase_map.csv` | significant site to PSP-annotated kinase, with coverage stated |
| `thesis/analysis/out/diffsites/diffsites_ksea.csv` | two-tailed KSEA on the differential result |
| `thesis/analysis/out/diffsites/diffsites_volcano.png` | log2FC against negative log10 BH-FDR, shortlist labelled |
| `thesis/analysis/out/diffsites/diffsites_sensitivity.csv` | the 12-per-group tier and the rank-test comparison |

The `.csv` and `.txt` outputs are small and would be committed. Estimated
runtime is a few minutes, dominated by reading the Excel file.

---

## 8. What I need before running

1. **Approval to run.** Explicitly requested and not assumed.
2. **A decision on §3's moderated-test question**, Welch now or wait for a
   moderated implementation.
3. **Confirmation of the shortlist rule in §5**, ideally checked with the
   collaborator, since they know what they can realistically assay.
4. Optional but valuable: **the patient pairing**, per §1. The analysis does not
   block on it, but the paired version is more powerful and the collaborator's
   confidence in a target would be higher.

---

## 9. Known limitations to state in any write-up

- **Provenance unresolved.** We still cannot cite a source study for this
  dataset. Gao et al. (2019) is ruled out on sample count, 159 patients against
  our 18 and 18.
- **One cohort, no external replication.** Nothing here is validated against a
  second dataset. For a journal submission, checking the top sites against a
  public HCC phosphoproteomics cohort would strengthen the result considerably,
  and is a separate piece of work.
- **34.7% missingness** means per-site power varies a lot. The observed counts
  are reported for every site so this stays visible.
- **1,319 multi-site rows are excluded**, matching the existing pipeline's
  behaviour. A site that only ever appears on a multiply-phosphorylated peptide
  is invisible to this analysis. If the collaborator cares about
  multiply-phosphorylated peptides, that is a deliberate extension, not a bug
  fix.
