# Differential phosphosites in HCC tumour vs adjacent non-tumour tissue

Candidate list for experimental validation. Prepared 2026-10-26.

Accompanying files: **`differential_shortlist.csv`** (659 sites, 15 columns) and
**`GRAPH_WALKTHROUGH.md`** (how the data and the network were built).

**Method.** Phosphosite abundances from 18 tumour and 18 adjacent non-tumour
samples, acquired in six batches of three and three. Values were log2
transformed and median centered per channel, because the supplied data are not
normalized: channel medians span a factor of 6.1. Each site was tested with a
linear model of log2 abundance on condition plus batch, necessary here because
the batch effect is large and changes direction between batches, with
Benjamini-Hochberg FDR across all tested sites. Nothing was imputed.

| Stage | Count |
|---|---|
| Phosphosites tested (>= 6 observations in both groups) | 2,973 |
| Significant at FDR < 0.05 | 1,548 |
| **Shortlist: FDR < 0.05, \|log2FC\| >= 1, >= 9 observations per group** | **659** |

Of the shortlist, 356 are down in tumour and 303 up, across 516 distinct genes.
550 (83.5%) stay significant without the batch term and 519 (78.8%) pass a
stricter 12-per-group filter. Both columns are in the CSV, and a site that
survives all of them is the safer bet for assay development.

## Top 10 by statistical confidence

All ten are down in tumour and all ten remain significant without the batch term.

| Site | log2FC | FDR | Obs. ctrl/tumour |
|---|---|---|---|
| PGM5-S122 | -3.59 | 2.3e-19 | 18 / 18 |
| PGLS-S49 | -2.59 | 2.6e-17 | 18 / 17 |
| HSPA12B-S25 | -3.52 | 1.3e-16 | 18 / 12 |
| ZPBP2-S173 | -1.90 | 7.8e-16 | 18 / 18 |
| PRKAR2B-S114 | -3.78 | 6.8e-15 | 18 / 18 |
| ASPDH-S168 | -2.76 | 1.4e-14 | 18 / 17 |
| DMTN-S372 | -1.92 | 1.7e-14 | 18 / 15 |
| TNKS1BP1-S1008 | -2.40 | 1.5e-13 | 18 / 12 |
| EML3-T881 | -2.87 | 2.5e-13 | 18 / 18 |
| ARHGAP31-S1178 | -1.26 | 4.9e-13 | 18 / 18 |

The largest changes by magnitude rather than confidence are CYP2A6-S131
(-5.52), PGM1-S509 (-5.36) and SQSTM1-S176 (+4.91), useful if the assay needs a
large signal.

## Kinase mapping

Only **111 of the 659 shortlist sites (16.8%)** have a kinase annotated to that
exact residue in PhosphoSitePlus. The other 83% are real measurements with no
annotated kinase, which reflects how incomplete the annotation is rather than
anything about the sites, and they are kept in the CSV. The kinases annotated to
the most significant sites are PRKACA (37), CDK1 (29), CSNK2A1 (22), MAPK1 (19)
and MAPK3 (18). Independently, activity inference across all measured substrates
(KSEA) finds 15 of 64 scoreable kinases significant, most decreased being
PRKACA, RPS6KA1, PRKD1, PAK1 and PRKG1. PRKACA tops both routes.

## Limitations, please read before selecting targets

**1. Not corrected for protein abundance.** A phosphosite can change because the
protein carrying it changed, not because its phosphorylation state changed. We
cannot correct for this: the protein expression sheet in the same file has 444
rows, is 91.7% empty, and **0 of the 659 shortlist sites** have usable protein
measurements. Any site taken forward should have its protein level checked
independently, for example by western blot alongside the phospho-specific assay.

**2. Unpaired test on a probably paired design.** The channel layout fits 18
matched pairs but the file does not record which control belongs to which
patient. This costs power rather than creating false positives, and we would
repeat the analysis as paired if the sample-to-patient table can be obtained.

**3. One cohort, 18 against 18, with no external replication, and the source
study is unidentified**, so we cannot retrieve its clinical annotation.

## Suggested filters

A conservative selection from the CSV: `BH_FDR` < 0.01, `n_control` >= 12 and
`n_tumour` >= 12, `also_sig_welch` true, and `n_psp_kinases` >= 1 if validating
the kinase as well as the substrate. Dropping that last filter gives many more
substrates but no kinase hypothesis to test alongside them.

Full output, including the volcano plot and every tested site, is in
`thesis/analysis/out/diffsites/`. To reproduce: `check_normalization`,
`differential_sites` and `diffsites_limitations` under
`.venv312/Scripts/python -m thesis.analysis.`
