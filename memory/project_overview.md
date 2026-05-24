---
name: project-overview
description: High-level architecture of the Kinase-Substrate HGN project — two pipelines, their outputs, and the biological goal
metadata:
  type: project
---

This is a master's-level research project building a heterogeneous graph network (HGN) for kinase-substrate phosphorylation prediction in a liver cancer context.

**Two pipelines:**

1. `src/` — Generic graph construction
   - Builds HGN from: Kinase_Substrate_Dataset, PPase substrates, PTMsigDB, MsigDB, PPI (high_confidence_score), PTMcode2
   - Entry: `python -m src.main`
   - Outputs: `outputs/nodes.csv.gz`, `outputs/edges.csv.gz`
   - Node types: protein (kinase/phosphatase/both/protein), site
   - Edge types: has_site, phosphorylates, dephosphorylates, site_same_pathway, protein_same_pathway, ppi_high_confidence, site_coevolution

2. `src/liver_network.py` — Liver extension
   - Expands generic graph with liver proteomics (LiverCancer_ProtExp_Phospho_casecntrl.xlsx)
   - Adds liver-specific nodes, reconnects via all databases, adds fold-change correlation edges
   - Entry: `python -m src.run_liver`
   - Outputs: `outputs/Liver_network_nodes.csv.gz`, `outputs/Liver_network_edges.csv.gz`
   - Fold-change correlation: healthy sample mean → tumor fold change → pairwise Pearson → 80th percentile threshold
   - Protein correlation disabled (only 3 usable protein rows in data — awaiting professor clarification)

3. `src_prediction/` — LOO-CV prediction pipeline
   - Evaluation set prep: `run_prepare_eval_set.py` → candidate kinases + liver positive edges
   - Freeze trials: `run_freeze_trials.py` → frozen_trials_10.csv.gz (debug: 10 trials)
   - Run LOO: `run_baseline_similarity.py` → node2vec embeddings + logistic regression scorer
   - Compare: `compare_experiments.py` → generic vs liver comparison

**Biological goal:** Identify liver-cancer-relevant kinase activity from phosphoproteomics. Validate that liver-specific graph context improves kinase ranking over generic graph. Ultimately discuss frequently predicted kinases and their liver cancer relevance; future work = experimental validation.

**Key scale numbers:**
- Generic graph: ~13,318 nodes, ~1,095,542 edges
- Liver graph: ~16,605 nodes, ~2,110,777 edges
- Candidate kinases: ~500 (from generic graph)
- Positive LOO edges (liver sites): see prep_summary.txt
- TTK has 68 substrates in the graph
