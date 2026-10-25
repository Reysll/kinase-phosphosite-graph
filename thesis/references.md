# References — verification status

Every entry below was checked against a publisher or PubMed Central record.
`VERIFIED` means the DOI resolves and the bibliographic details were confirmed.
`UNVERIFIED` entries must not be cited until resolved.

APA 7 formatting will be applied at the docx build; the list here is ordered by
topic for checking, not alphabetically.

## Annotation gap and the dark kinome

**VERIFIED** — Deznabi, I., Arabaci, B., Koyutürk, M., & Tastan, O. (2020).
DeepKinZero: Zero-shot learning for predicting kinase–phosphosite associations
involving understudied kinases. *Bioinformatics, 36*(12), 3652–3661.
https://doi.org/10.1093/bioinformatics/btaa013 · PMC7320620
*Supports:* ">95% of reported human phosphosites have no known kinases"
(verbatim quote confirmed).

**VERIFIED** — Maier, B. D., Petursson, B., Lussana, A., & Petsalaki, E. (2025).
Data-driven extraction of human kinase-substrate relationships from omics
datasets. *Molecular & Cellular Proteomics, 24*(8), 100994.
https://doi.org/10.1016/j.mcpro.2025.100994 · PMC12359690
*Supports:* "~150 kinases lack known substrates, and ~80% of annotated
phosphosites are attributed to just 20% of the kinases" (verbatim).

**VERIFIED** — Needham, E. J., Parker, B. L., Burykin, T., James, D. E., &
Humphrey, S. J. (2019). Illuminating the dark phosphoproteome.
*Science Signaling, 12*(565), eaau8645.
https://doi.org/10.1126/scisignal.aau8645
*Note:* original source of the 80/20 statistic that Maier et al. cite.

## Data resources

**VERIFIED** — Hornbeck, P. V., Kornhauser, J. M., Tkachev, S., Zhang, B.,
Skrzypek, E., Murray, B., Latham, V., & Sullivan, M. (2012). PhosphoSitePlus: A
comprehensive resource for investigating the structure and function of
experimentally determined post-translational modifications in man and mouse.
*Nucleic Acids Research, 40*(D1), D261–D270.
https://doi.org/10.1093/nar/gkr1122

**VERIFIED** — Duan, G., Li, X., & Köhn, M. (2015). The human DEPhOsphorylation
database DEPOD: A 2015 update. *Nucleic Acids Research, 43*(D1), D531–D535.
https://doi.org/10.1093/nar/gku1009

**VERIFIED** — Krug, K., Mertins, P., Zhang, B., Hornbeck, P., Raju, R., Ahmad,
R., … Carr, S. A. (2019). A curated resource for phosphosite-specific signature
analysis. *Molecular & Cellular Proteomics, 18*(3), 576–593.
https://doi.org/10.1074/mcp.TIR118.000943
*Note:* the PTMsigDB resource paper.

**VERIFIED** — Liberzon, A., Birger, C., Thorvaldsdóttir, H., Ghandi, M.,
Mesirov, J. P., & Tamayo, P. (2015). The Molecular Signatures Database (MSigDB)
hallmark gene set collection. *Cell Systems, 1*(6), 417–425.
https://doi.org/10.1016/j.cels.2015.12.004

**VERIFIED** — Minguez, P., Letunic, I., Parca, L., Garcia-Alonso, L., Dopazo,
J., Huerta-Cepas, J., & Bork, P. (2015). PTMcode v2: A resource for functional
associations of post-translational modifications within and between proteins.
*Nucleic Acids Research, 43*(D1), D494–D502.
https://doi.org/10.1093/nar/gku1081

**UNVERIFIED — DO NOT CITE YET** — the protein-protein interaction source.
The repository supplies `data/raw/high_confidence_score.csv` with no provenance
record. It is commonly assumed to be STRING, but nothing in the code or docs
states this. Either establish the origin or describe it generically as
"a high-confidence protein interaction set" with the file documented in an
appendix.

## Graph representation learning

**VERIFIED** — Grover, A., & Leskovec, J. (2016). node2vec: Scalable feature
learning for networks. *Proceedings of the 22nd ACM SIGKDD International
Conference on Knowledge Discovery and Data Mining*, 855–864.
https://doi.org/10.1145/2939672.2939754

**VERIFIED** — Perozzi, B., Al-Rfou, R., & Skiena, S. (2014). DeepWalk: Online
learning of social representations. *Proceedings of the 20th ACM SIGKDD
International Conference on Knowledge Discovery and Data Mining*, 701–710.
https://doi.org/10.1145/2623330.2623732

**VERIFIED** — von Luxburg, U. (2007). A tutorial on spectral clustering.
*Statistics and Computing, 17*(4), 395–416.
https://doi.org/10.1007/s11222-007-9033-z

**VERIFIED** — Shi, J., & Malik, J. (2000). Normalized cuts and image
segmentation. *IEEE Transactions on Pattern Analysis and Machine Intelligence,
22*(8), 888–905. https://doi.org/10.1109/34.868688

**PARTIAL** — Ng, A. Y., Jordan, M. I., & Weiss, Y. (2002). On spectral
clustering: Analysis and an algorithm. *Advances in Neural Information
Processing Systems, 14*, 849–856.
*Note:* NeurIPS proceedings of this era have no DOI. Freely accessible at
https://proceedings.neurips.cc/paper/2001/hash/801272ee79cfde7fa5960571fee36b9b-Abstract.html
Acceptable under the "freely accessible link" rule, but flag for the committee.

## Kinase activity inference

**VERIFIED** — Wiredja, D. D., Koyutürk, M., & Chance, M. R. (2017). The KSEA
App: A web-based tool for kinase activity inference from quantitative
phosphoproteomics. *Bioinformatics, 33*(21), 3489–3491.
https://doi.org/10.1093/bioinformatics/btx415 · PMC5860163

## Liver cancer dataset

**NOT CITED — source study unidentified.** The liver dataset
`data/liver/LiverCancer_ProtExp_Phospho_casecntrl.xlsx` is described by its
contents and sample counts throughout the thesis, with no attribution. As of
2026-10-25 it is cited nowhere in Chapters II or III.

Evidence gathered so far, for whoever resolves this:

- The sheets retain their original captions. `ProteinExpression` is "Table S4
  Detailed information of unique proteins identified in HCC tumour tissues" and
  `Phosphorylation` is "Table S11 Identified and quantified phosphopeptides,
  phosphosites and corresponding phosphoproteins in HCC tumour and adjacent
  non-tumour tissues". So the source is some HCC proteogenomics paper whose
  supplementary tables S4 and S11 these are.
- British spelling ("tumour") suggests a non-US journal.
- The file holds **18 tumor + 18 control** samples in six batches of three plus
  three, with 444 protein rows and 6,571 phosphopeptide rows.
- **Gao et al. (2019), *Cell* 179(2), 561–577.e22,
  https://doi.org/10.1016/j.cell.2019.08.052, is RULED OUT**: that study reports
  a cohort of 159 patients. The publication itself verifies cleanly; it is
  simply not the source of this file.
- Two extended literature searches on the captions did not identify the paper.
  The next step is to ask the advisor or the file's supplier, not to search
  again.

## Hepatocellular carcinoma biology (Chapter VI, held)

**VERIFIED** — Jiang, X., et al. (2023). A RIPK3-independent role of MLKL in
suppressing parthanatos promotes immune evasion in hepatocellular carcinoma.
*Cell Discovery, 9*, 7. https://doi.org/10.1038/s41421-022-00504-0 · PMC9845215

**VERIFIED** — Köles, L., Ribiczey, P., Szebeni, A., et al. (2024). The role of
TRPM7 in oncogenesis. *International Journal of Molecular Sciences, 25*(2), 719.
https://doi.org/10.3390/ijms25020719 · PMC10815510

**PARTIAL** — Necroptosis-driving genes RIPK1, RIPK3 and MLKL-p are associated
with intratumoral CD3+ and CD8+ T cell density and predict prognosis in
hepatocellular carcinoma. PMC8915343. Full author list, journal and DOI still
to be captured.

## Summary

| Status | Count |
|---|---|
| VERIFIED with DOI | 14 |
| PARTIAL (accessible link, missing DOI or details) | 2 |
| UNVERIFIED, not cited | 1 (PPI source) |
| Source unidentified, not cited | 1 (liver dataset) |

## In-text citation format (APA 7)

Applied to Chapter II on 2026-10-25:

- Parenthetical citations use `&`: (Grover & Leskovec, 2016).
- Narrative citations use "and" in prose and place the year in parentheses:
  Deznabi et al. (2020) state that...
- Three or more authors are always `et al.` from the first citation, so
  "Ng, Jordan, & Weiss, 2002" became "Ng et al., 2002".
- Multiple works in one parenthesis are ordered alphabetically and separated by
  semicolons: (Krug et al., 2019; Liberzon et al., 2015).
- Working `[VERIFIED]` and `[PARTIAL]` markers have been removed from the prose.
  This file is the verification record.
