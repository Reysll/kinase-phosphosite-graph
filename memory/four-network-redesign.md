---
name: four-network-redesign
description: Full project history — embedding leakage fix, four networks, spectral embedding (Phase 2), category-based ranking, current results
metadata:
  type: project
---

## Phase 1 — node2vec, 6,909 trials (complete 2026-05-27)

- Leakage fix: PSP `phosphorylates` edges stripped from embedding graph.
- `adjusted_held_out_rank`: primary metric (co-kinases removed before ranking).
- Four networks: generic, control (healthy), cancer (tumor), liver FC.

Results:
| Network | adj median rank |
|---|---|
| generic | 188.0 |
| control | 166.5 |
| cancer  | 174.0 |
| liver FC| 165.0 |

Hub bias: ULK1=100% top-1 in Liver FC; VRK2=94.8% in Control; VRK3=95.5% in Cancer.

## Phase 2 — SpectralEmbeddingStrategy + category ranking (complete 2026-06-17)

### Option 1: SpectralEmbeddingStrategy
- File: `src_prediction/embedding_strategy.py`
- Normalized Laplacian eigenvectors. directed=False. scipy ARPACK eigsh, 32 components.
- Implemented instead of GCL (Dr. Ayati's suggestion) — same degree-normalization, closed-form.

### Option 2: Category-based ranking
- File: `src_prediction/kinase_categories.py`
- Thresholds: poor ≤5 (166 kinases), average 5-20 (133), rich >20 (121).
- Pre-filter PSP substrate counts from full generic graph.
- Hub placement: VRK3→poor, VRK2→average, ULK1→rich.
- New columns: `held_out_kinase_category`, `adjusted_held_out_rank_in_category`, `top1_poor/average/rich`.

### Full results (6,909 trials, spectral + categories)

|                  | Generic | Control | Cancer  | Liver FC |
|---|---|---|---|---|
| adj mean rank    | 126.6   | 105.1   | 105.0   | 105.1    |
| adj median rank  | 79.0    | 74.0    | 74.0    | 74.0     |
| adj top-1        | 0.2%    | 2.7%    | 0.1%    | 0.0%     |
| adj top-20       | 10.6%   | 22.1%   | 22.1%   | 22.0%    |
| best-true median | 34.0    | 15.0    | 15.0    | 15.0     |

vs Phase 1 node2vec adj median: Generic 188→79 (-58%), Liver FC 165→74 (-55%).

### Key findings
1. Spectral halves median rank across all networks.
2. Control=Cancer=Liver FC (all median 74) under spectral — correlation TYPE no longer differentiates. With node2vec, they differed significantly.
3. Generic clearly weakest (79 vs 74) — correlation edges add structural signal.
4. Control has highest top-1 (2.7% vs ~0%). Most distributed global top-1: CDK1 57% / PIK3R1 43%.
5. Hub kinases changed: from VRK family to CDK1, PIK3R1, TGFBR2 — more biologically plausible.
6. Within-category hub bias persists (100% concentration per category). Categories reduce cross-category competition but not embedding-level bias.

### Top-1 global monopoly (spectral)
| Network  | Top-1 hub |
|---|---|
| Generic  | TTK 92.9%, IRAK4 7.1% |
| Control  | CDK1 57.1%, PIK3R1 42.9% |
| Cancer   | PIK3R1 99.3% |
| Liver FC | PIK3R1 100.0% |

### Per-category liver spectral (6,909 trials)
| Category | n trials | adj median rank | within-cat top-1 | hub |
|---|---|---|---|---|
| poor     | 201      | 79              | 1.0%  | PIK3R1 100% |
| average  | 836      | 61              | 0.4%  | TGFBR2 100% |
| rich     | 5872     | 32              | 4.8%  | CDK1 100%   |

## Phase 3 — analyze_results.py + KSEA (complete 2026-06-17)

### analyze_results.py extended (8 analyses)
- Analysis 6: `phase_comparison.txt` — node2vec vs spectral for all 4 networks
- Analysis 7: `category_ranks.txt` — per-category (poor/average/rich) rank breakdown from spectral
- Analysis 8: `ranking_shifts_spectral.txt` — generic vs liver FC delta under spectral
- TTK deepdive now covers both Phase 1 (node2vec) and Phase 2 (spectral)
- New constants: `SPECTRAL_EXPERIMENT_DIRS`, `SPECTRAL_EXPERIMENT_LABELS`, `PHASE_PAIRS`

### KSEA (src_prediction/ksea.py + run_ksea.py)
- Formula: Wiredja et al. 2017 (Bioinformatics 33:3489) — score = (s_bar - p_bar) / (m * delta)
  where m is substrate count (in denominator — penalizes kinases with many substrates).
  P-values: one-tailed N(0,1) + Benjamini-Hochberg FDR (scipy.stats.false_discovery_control).
- FC from Phosphorylation sheet: log2(mean_tumor / mean_ctrl), 3,228 sites
- PSP KSEA: 64 kinases scored (≥3 substrates in FC dataset). No significant hits (adj_p ~0.5).
  - Top-ranked: MAPK12 (n=3, score=1.08), CDKL5 (n=3, score=1.05) — few substrates, very high FC
  - CDK1 (n=61, score=0.004), CSNK2A1 (n=35, score=0.015) — penalized by large m
  - Dataset underpowered: m in denominator means large-substrate kinases can't score well
- Predicted substrate KSEA: hub-dominated — TTK (generic, ~150 sites) and PIK3R1 (liver, ~164 sites)
  Both score negative (their "predicted substrates" are random mix — no FC enrichment).
- Key limitation: LOO results only give top-1 per site; full-inference pass needed for richer substrate sets.

## Next steps
- Draft email to Dr. Ayati (Phase 2 results + KSEA findings).
- Consider "predict all" inference pass (score all kinases for all liver FC sites) for richer KSEA.
- Biological interpretation: CDK1/CSNK2A1 as top KSEA hits in liver cancer context.
