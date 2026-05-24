---
name: project-status
description: Current state of the project as of 2026-05-21 — what works, results so far, what's blocked, next steps
metadata:
  type: project
---

**Date:** 2026-05-21. Project at debug/proof-of-concept stage (10 frozen trials).

**Already resolved (per email thread with professor):**
- Graph expansion: liver nodes added and reconnected via all databases ✓
- Multi-true-kinase evaluation: reports both held-out rank AND best-true rank ✓
- Scorer: upgraded from cosine similarity to logistic regression on node2vec pair features ✓
- Site fold-change correlation: healthy mean → tumor FC → Pearson → 80th percentile threshold ✓

**Current 10-trial results (generic vs liver, 80% threshold):**
- Generic: mean held-out rank 146.4, median 118.0
- Liver 80%: mean held-out rank 73.3, median 58.0 → ~2× improvement, 7/10 trials improved
- Liver 85% was worse (mean 144.8, median 128.5) — 80% is the working setting

**TTK analysis (from debug run):**
- TTK has 68 substrates in the graph
- Liver model predicts TTK as top-1 in 3/10 debug trials (generic: 1/10)
- Professor wanted more TTK hits to argue TTK importance in liver cancer — not achieved at debug scale

**Performance fix completed (2026-05-21):**
- `leave_one_out.py`: node2vec now runs ONCE before trial loop (was 50× per trial)
- `pair_features.py`: `build_pair_feature_table` vectorized with numpy matrix ops (was itertuples loop)
- `TrialTask` now carries pre-built `embeddings` dict instead of `graph_edges`/`node2vec_params`
- Public API of `run_leave_one_out` unchanged — callers need no changes

**Open/blocked:**
- Protein co-expression: only 3 usable protein rows → disabled; professor has not clarified data source
- Full LOO run not yet completed (was blocked on performance — now unblocked)

**Next steps:**
1. Run full LOO on liver graph (all available trials, not just 10)
2. TTK deep-dive: for all 68 TTK substrates, report top-k accuracy in full run
3. Sites with multiple known kinases: identify and check biology (professor's explicit ask)
4. Most frequent top-1 predicted kinases in liver model vs generic
5. Extended abstract: liver-specific context → improved predictions → liver cancer kinase biology → future experimental validation
