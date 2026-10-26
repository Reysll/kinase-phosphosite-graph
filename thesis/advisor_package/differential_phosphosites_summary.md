# Differential phosphosites in HCC tumour vs adjacent non-tumour tissue

Candidate list for experimental validation. Prepared 2026-10-26.

Accompanying file: **`differential_shortlist.csv`**, 659 sites, 15 columns.

## What was done

Phosphosite abundances from 18 tumour and 18 adjacent non-tumour samples,
acquired in six batches of three tumour and three control each. Abundances were
log2 transformed and median centered per channel, because the supplied values
are not normalized: channel medians spanned a factor of 6.1.

Each site was tested with a linear model of log2 abundance on condition plus
batch, which is necessary here because the batch effect is large and changes
direction between batches. Benjamini-Hochberg FDR across all tested sites.
No values were imputed; each site is tested on the channels where it was
actually observed.

## Numbers

| Stage | Count |
|---|---|
| Phosphosites tested (>= 6 observations in both groups) | 2,973 |
| Significant at FDR < 0.05 | 1,548 |
| **Shortlist: FDR < 0.05, \|log2FC\| >= 1, >= 9 observations per group** | **659** |

Within the shortlist:

| Property | Count |
|---|---|
| Down in tumour | 356 |
| Up in tumour | 303 |
| Distinct genes | 516 |
| Also significant without the batch term | 550 (83.5%) |
| Also pass the stricter >= 12 per group filter | 519 (78.8%) |

The last two rows are robustness checks. A site that survives all of them is a
safer bet for assay development.

## Top 10 by statistical confidence

Ordered by FDR. All ten are down in tumour and all ten remain significant
without the batch term.

| Site | log2FC | FDR | Obs. (ctrl/tumour) |
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

Three of these carry a known kinase annotation: PRKAR2B-S114 (PRKACA),
EML3-T881 (CDK1) and ARHGAP31-S1178 (RPS6KA1).

## Largest changes

Ordered by effect magnitude rather than confidence. Useful if the assay needs a
large signal.

| Site | log2FC | FDR | Annotated kinase |
|---|---|---|---|
| CYP2A6-S131 | -5.52 | 4.3e-12 | none |
| PGM1-S509 | -5.36 | 2.3e-07 | none |
| SQSTM1-S176 | **+4.91** | 1.5e-03 | none |
| HSPB6-S16 | -4.30 | 7.9e-10 | PRKACA, PRKD1 |
| PRKAR2B-S114 | -3.78 | 6.8e-15 | PRKACA |

The direction is coherent with liver biology: CYP2A6 and PGM1 are
hepatocyte-function proteins losing phosphorylation in tumour, while SQSTM1
(p62) gains it.

## Kinase mapping

Two routes were used, and the coverage of the first is the main caveat.

**Direct annotation.** Only **111 of the 659 shortlist sites (16.8%)** have any
kinase annotated to that exact residue in PhosphoSitePlus. The other 83% are
real measurements with no annotated kinase, which reflects how incomplete the
annotation is rather than anything about the sites. They are kept in the CSV.
The kinases most often annotated to a significant site are PRKACA (37 sites),
CDK1 (29), CSNK2A1 (22), MAPK1 (19) and MAPK3 (18).

**Activity inference (KSEA).** Run on the fold changes across all measured
substrates, 15 of 64 scoreable kinases reach FDR < 0.05. Most decreased:
PRKACA, RPS6KA1, PRKD1, PAK1, PRKG1. Most increased: CDK2, CSNK2A1, CDKL5,
CDK7, CDK1.

PRKACA appears at the top of both routes, which is the most consistent kinase
signal in the data.

## Limitations, please read before selecting targets

**1. These are not corrected for protein abundance.** A phosphosite can change
because the protein carrying it changed, not because its phosphorylation state
changed. We cannot correct for this with the data we have: the protein
expression sheet in the same file has 444 rows, is 91.7% empty, and **0 of the
659 shortlist sites** have protein measurements sufficient to support a
correction. Any site taken forward should have its protein level checked
independently, for example by western blot alongside the phospho-specific assay.

**2. The test is unpaired although the design is probably paired.** The channel
layout is consistent with 18 matched tumour and adjacent-normal pairs, but the
file does not record which control channel belongs to which patient. An unpaired
test is conservative, so this costs power rather than creating false positives.
If the sample-to-patient table can be obtained, the analysis would be repeated
as paired.

**3. One cohort, no external replication.** 18 against 18, with no validation in
a second dataset.

**4. The source study is unidentified.** The file's sheets are supplementary
tables S4 and S11 of some HCC proteogenomics paper that we have not been able to
identify, so we cannot cite it or retrieve its clinical annotation.

## Suggested filters for the collaborator

Columns are in the CSV. A conservative selection would be:

- `BH_FDR` < 0.01
- `n_control` >= 12 and `n_tumour` >= 12
- `also_sig_welch` is true
- `n_psp_kinases` >= 1, if validating the kinase as well as the substrate

Relaxing the last filter gives many more substrates but no kinase hypothesis to
test alongside them.

## Reproducing this

All read-only, writing nothing outside `thesis/`.

```
.venv312/Scripts/python -m thesis.analysis.check_normalization
.venv312/Scripts/python -m thesis.analysis.differential_sites
.venv312/Scripts/python -m thesis.analysis.diffsites_limitations
```

Full output, including the volcano plot and every tested site rather than just
the shortlist, is in `thesis/analysis/out/diffsites/`.
