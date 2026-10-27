# Why the network context never reaches the prediction

Written 2026-10-26. Read-only diagnosis of existing outputs and embeddings.
No project code, data or results were modified.

This document explains a single defect that accounts for most of the negative
results in this project, and separates what has been measured from what is
still hypothesis. It supersedes the K-block explanation offered in
`CATCH_UP.md` (c), which was a reasonable guess and is not what is happening.

> **Update 2026-10-27: the fix was built, run, and it did not work.** The
> diagnosis below is confirmed on the embedding and **refuted as an explanation
> of the ranking**. Discarding the whole null space removes the degree
> dimensions (17 of 32 down to 0 of 32) and leaves the predicted ranking
> degree-ordered at Spearman rho = -0.9797. The median held-out rank does not
> move. Section 9 has the numbers and the measured reason. Read Section 9
> before acting on Sections 3 or 7.

---

## 1. The short version

The model is not using the graph context. It is sorting kinases by their degree
in the graph, descending, and returning that same order at every phosphosite.

The correlation between a kinase's embedding-graph degree and its mean
predicted rank is **Spearman rho = -0.9995** under Control and **-0.9973**
under Liver FC. For comparison, PSP substrate count, which is the popularity
baseline's ranking key, reaches only -0.44.

There is a bug in how the spectral embedding is computed.
`SpectralEmbeddingStrategy` discards one eigenvector where it should discard
the whole null space of the Laplacian, and the vectors it fails to discard are
proportional to sqrt(degree). Between **12 and 17 of the 32 embedding
dimensions** are therefore rescaled copies of node degree.

**That bug is real but it is not the cause of the ranking.** Removing those
dimensions leaves the ranking degree-ordered (rho = -0.9797, from -0.9995).
The degree signal does not live only in the null space: every eigenvector of
the normalized Laplacian is `u = D^1/2 v`, so each node's embedding has norm
proportional to sqrt(degree), and that scale survives any choice of which
eigenvectors to keep. Measured kinase norm against degree: rho = +0.94 both
before and after the fix. Section 9.

---

## 2. What the model is actually doing

From the Phase 3 predict-all output already on disk, for each kinase take its
mean `rank_in_site` across all scored sites, then regress that on the kinase's
degree in the embedding graph.

| Network | Embedding | Degree rho | Degree R^2 | Substrates rho | Rank SD |
|---|---|---|---|---|---|
| Generic | **spectral** | **-0.9387** | 0.8546 | -0.3272 | 4.93 |
| Generic | node2vec | -0.3206 | 0.1070 | -0.0285 | 10.20 |
| Generic | concat | -0.3252 | 0.1111 | -0.0500 | 9.71 |
| Control | **spectral** | **-0.9995** | 0.9707 | -0.4382 | **1.05** |
| Control | node2vec | -0.4362 | 0.1949 | -0.1227 | 9.08 |
| Control | concat | -0.4248 | 0.1833 | -0.1243 | 9.06 |
| Liver FC | **spectral** | **-0.9973** | 0.9669 | -0.4483 | **1.18** |
| Liver FC | node2vec | -0.3594 | 0.1354 | -0.1788 | 8.20 |
| Liver FC | concat | -0.3824 | 0.1466 | -0.1622 | 7.57 |

Degree is computed on the *embedding* graph, which has all `phosphorylates`
edges stripped, so it carries no direct information about the prediction
target. This is not label leakage. It is the embedding handing the classifier a
degree channel and the classifier using it.

**The defect is specific to the spectral embedding.** Node2vec and concat reach
only -0.32 to -0.44, against spectral's -0.94 to -1.00. That was not the
expected result: node2vec is generally held to encode degree through
random-walk frequency, and the prediction before measuring was that it would
look similar. It does not, which strengthens the attribution to the null-space
bug in Section 3 rather than to graph embedding in general.

The `Rank SD` column is the median per-kinase standard deviation of
`rank_in_site` across sites, and it points the other way from every ranking
metric this project has reported. Node2vec varies **8 to 10 rank positions**
per kinase across sites; spectral varies **1**. Spectral produces better median
ranks while being almost completely site-blind, because a degree ordering is a
serviceable global prior. Node2vec ranks worse and is the only one of the three
with ranks that move across sites at all.

**Caveat: part of node2vec's rank SD may be seed noise, so "node2vec is more
site-aware" is not established.** Every number in the table above comes from a
single seed. Findings 11 and 12 measured seed-to-seed variation under node2vec
that exceeds the context effect (+0.373 against +0.249), with six runs producing
six different top-1 hubs and a Generic median spread of 40 ranks. Node2vec's
random walks are resampled on every run, so an 8 to 10 position spread in a
kinase's rank across sites is consistent with site information, with run-to-run
walk noise, or with any mixture of the two. Separating them needs the same
kinase and site scored under several seeds, which has not been run. Spectral is
deterministic given the graph, so its near-constant ranks are not subject to this
caveat; what stands is that spectral's ranks barely move across sites, not that
node2vec's movement is biological.

The top and bottom of the Liver FC ranking make the point without statistics:

| Kinase | Mean rank | Rank SD | Degree | PSP substrates |
|---|---|---|---|---|
| PIK3R1 | 1.23 | 0.43 | 1,801 | **2** |
| EGFR | 2.12 | 0.82 | 1,787 | 101 |
| CDK1 | 2.74 | 0.64 | 1,789 | 693 |
| PRKCA | 4.12 | 0.42 | 1,705 | 502 |
| SGK1 | 4.81 | 0.39 | 1,688 | 54 |
| ... | | | | |
| STK26 | 413.0 | **0.00** | 10 | 13 |
| TSSK4 | 414.0 | **0.00** | 1 | 1 |
| MAP3K20 | 415.0 | **0.00** | 12 | 7 |

PIK3R1 is ranked first at essentially every one of 3,228 sites while having two
annotated substrates. It is the highest-degree node in the graph. CDK1, with
693 substrates, ranks third. The bottom entries have a rank standard deviation
of exactly zero, meaning they hold the identical rank at every site.

This also identifies the Phase 2 "hub" winners for what they were. PIK3R1 took
top-1 at 100% of Liver FC trials because it is the highest-degree node, not
because of anything about PI3K signalling in liver cancer.

---

## 3. Why: the Laplacian null space

`SpectralEmbeddingStrategy.fit` builds the normalized Laplacian
`L = I - D^-1/2 A D^-1/2` and then does:

```python
k = min(self.dimensions + 1, n - 1)
eigenvalues, eigenvectors = spla.eigsh(L, k=k, which="SM", v0=...)
order = np.argsort(eigenvalues)[1 : self.dimensions + 1]
```

It requests 33 eigenvectors and drops exactly one, on the assumption that the
graph is connected and so has a single eigenvalue at zero. **The embedding
graphs are not connected.**

| Network | Nodes | Components (size >= 2) | Giant component | Zero eigenvalues among the 33 |
|---|---|---|---|---|
| Generic | 13,313 | 34 | 13,223 | 13 |
| Control | 17,731 | 116 | 17,473 | 18 |
| Liver FC | 17,731 | 122 | 17,463 | 17 |

For a graph with `c` components the normalized Laplacian has eigenvalue 0 with
multiplicity `c`. The spectra confirm it: the numerically-zero eigenvalues sit
at 1.2e-08 to 1.8e-08 and the first genuinely non-zero one is 3.53e-02
(Generic), 3.65e-03 (Control) and 3.68e-03 (Liver FC), a gap of five orders of
magnitude, so the count is unambiguous.

Null space of `L` is spanned by the component indicators scaled by sqrt(degree):
for component `C`, the vector `u_C(i) = sqrt(deg_i)` on `C` and 0 elsewhere
satisfies `L u_C = 0`. ARPACK returns an arbitrary orthonormal basis of that
space, so restricted to the nodes of any single component, every basis vector
is a scalar multiple of sqrt(deg).

Measured directly, as R^2 against sqrt(deg) over giant-component nodes:

| Network | Shipped dims that are sqrt(deg) | Median R^2, shipped | After dropping the whole null space |
|---|---|---|---|
| Generic | **12 / 32** | 0.0012 | **0 / 32** |
| Control | **17 / 32** | **1.0000** | **0 / 32** |
| Liver FC | **16 / 32** | 0.5219 | **0 / 32** |

Generic's null dimension reaches R^2 = 1.0000 while its informative dimensions
sit at 0.0011, which is the cleanest possible separation. For Control the
*median* of all 32 shipped dimensions is R^2 = 1.0000, meaning at least half
that embedding is degree.

The fix is cheap: request enough eigenvectors and keep only those with
eigenvalue above a tolerance. Computing 72 eigenvectors took 0.1 minutes per
network, so this is not a performance tradeoff.

It is also not sufficient. Everything in this section is a correct statement
about the *embedding*; Section 9 shows that acting on it does not change the
*ranking*, because the degree signal is not confined to the null space.

A note on the docstring. `SpectralEmbeddingStrategy` claims "High-degree nodes
are penalized by construction" and that spectral embeddings do not share
node2vec's hub bias, and `CLAUDE.md` repeats the second claim. The
normalization argument is correct for the non-trivial eigenvectors. It does not
hold for the null space the code fails to discard, which is where the degree
information enters. Both claims are not merely unsupported but inverted:
spectral is the *more* degree-dependent of the two embeddings, at rho = -0.9995
against node2vec's -0.4362.

### 3.1 What degree is made of

Degree is not a neutral measure of biological connectivity here. Counting edges
that touch a candidate kinase, with `phosphorylates` excluded:

| Relation | Edges touching a kinase | Share |
|---|---|---|
| `protein_same_pathway` | 273,279 | **95.5%** |
| `ppi_high_confidence` | 11,115 | 3.9% |
| `has_site` | 1,824 | 0.6% |
| `dephosphorylates` | 1 | 0.0% |

So a kinase's degree is, to 95.5%, its MSigDB pathway-clique membership
weighted by clique size. That closes the chain back to the preprocessing audit:
`protein_same_pathway` is 90.6% of generic edges
(`GRAPH_WALKTHROUGH.md` section 3.3), those edges are an annotation rather than
a measurement, and they are what the ranking is built from.

The chain end to end:

    MSigDB pathway cliques (95.5% of kinase edges)
      -> kinase degree
      -> 12-17 of 32 spectral dimensions (Laplacian null space)
      -> predicted rank (rho = -0.9995)
      -> site-invariance, null-equivalence, hub bias, baseline dominance

---

## 4. The second defect: the sites are indistinguishable

Removing the degree channel is necessary but not sufficient. There is a second,
structural problem, and it is the one that actually blocks context-aware
prediction.

Measure the share of site-embedding magnitude that is common to every site:

    rho = ||S_bar|| / sqrt(||S_bar||^2 + E||S - S_bar||^2)

| Network | Site rho, as shipped | After dropping null space | After restricting to giant component |
|---|---|---|---|
| Generic | 0.1190 | 0.1239 | 0.1243 |
| Control | 0.2952 | 0.2795 | 0.2824 |
| Liver FC | **0.9160** | **0.8882** | **0.8882** |

Under Liver FC, 89% of the site-embedding magnitude is identical across all
3,228 sites even after both fixes. The sites are nearly the same point. A
consequence, measured on the `K*S` block of the feature vector: 95% of its
magnitude is a per-kinase constant that no choice of logistic-regression
weights can make site-dependent (median over kinases; 0.88 for Control, 0.47
for Generic).

The spectral gap is where the two regimes differ. Generic's first non-zero
eigenvalue is 3.53e-02; the context networks' are 3.65e-03 and 3.68e-03, an
order of magnitude smaller. Adding roughly 976,855 correlation edges leaves the
graph with no bottleneck structure for an embedding to resolve.

---

### 4.1 Which correlation filter would fix it

Testable on the embedding alone, with no model run. Each variant keeps the same
base relations and a filtered subset of the Liver FC correlation edges, with
the whole null space discarded so differences are structural. Generic's
reference site rho is 0.1239.

| Filter | Corr edges | Components | Gap | Site rho | Kinase rho | Sites kept |
|---|---|---|---|---|---|---|
| Current, 80th pct | 976,855 | 122 | 3.68e-03 | 0.8882 | 0.6938 | 3,228 |
| FDR < 0.05 | 846,586 | 122 | 3.93e-03 | 0.8462 | 0.6824 | 3,228 |
| n_common >= 12 | 379,783 | 338 | 5.53e-03 | 0.2244 | 0.9189 | 3,228 |
| FDR and n >= 12 | 368,439 | 339 | 5.62e-03 | 0.2285 | 0.9184 | 3,228 |
| **Top 25% by \|r\|** | 244,214 | **127** | 8.08e-03 | 0.3685 | **0.5632** | **3,228** |
| n_common = 18 | 78,138 | 667 | 1.14e-02 | 0.1622 | 0.9433 | 2,529 |
| No correlation edges | 0 | 998 | 3.02e-02 | 0.0588 | 0.8141 | 1,860 |

Three things follow.

**The collapse is caused by the correlation edges and is reversible.** Site rho
falls monotonically as they are pruned, from 0.8882 to 0.0588.

**Significance testing is the wrong filter.** FDR < 0.05 removes only 13% of
the edges and moves site rho by 0.04, because the pool has already been cut at
the pooled 80th percentile, so 86.7% of what remains passes FDR. This matters
because adding FDR control is the natural-sounding fix and the measurement says
it would not work. Sample count and magnitude are what bind.

**A magnitude threshold is preferable to a sample-count floor**, which was not
the expected answer. A sample-count floor removes edges systematically: a site
detected in few samples loses *all* of its correlation edges and is stranded,
so components rise from 122 to 338 to 667 and site coverage starts falling. A
magnitude threshold prunes uniformly, holding components at 127 and keeping all
3,228 sites, while giving the largest gap of any full-coverage variant and the
best kinase rho of any variant at all.

This also quantifies the warning in `CLAUDE.md` about using base relations only
on the liver graph. It is accurate: 998 components, and 1,368 of 3,228 sites
(42%) stranded outside the giant component. The eigensolver "breakage" it
describes is the null-space multiplicity of Section 3, which the corrected
eigensolver handles, so the variant is now measurable even though it remains a
poor choice on coverage grounds.

Note the tradeoff in the kinase column. Every filter that de-collapses the
sites except the magnitude one pushes kinase rho up, to 0.92-0.94, trading site
degeneracy for kinase degeneracy. Any filter chosen has to be checked on both.

## 5. What this explains

Four findings previously treated as separate reduce to this one cause.

**Finding 9, the degree-preserving null reproduces 96% of ranks.** This stops
being a surprising empirical result and becomes close to an identity. A
degree-preserving null preserves the degree sequence exactly, and the ranking
is a function of degree, so the ranking is preserved. The null experiment was
correctly executed and correctly interpreted as "the biology is not being
used"; what it could not reveal was *why*.

**Findings 13 and 26, eta^2 = 0.9997 rising to 1.0000.** Degree is a per-kinase
constant, so a degree-ordered ranking is site-invariant by construction. This
also explains why fixing the negative sampler pushed eta^2 *up*: cleaner
negatives let the classifier fit the degree direction more precisely.

**Finding 14, the popularity baseline beats every model.** The baseline ranks by
substrate count, which is weakly predictive. The models rank by degree, which
is not predictive. These are different orderings, correlating at only
rho = -0.44, so the baseline is not a better version of the model. It is
measuring something real while the model measures connectivity.

**Phase 2's hub kinases.** PIK3R1, CDK1 and TGFBR2 were read as cancer-relevant
nodes replacing node2vec's VRK family. They are the highest-degree nodes.

**The Phase 2 choice of embedding was backwards.** Spectral was adopted over
node2vec because it halved the median rank and was believed to fix the hub
bias. On the evidence here it halved the median rank *by* encoding a degree
prior, and it is the more degree-dependent of the two by a wide margin
(rho = -0.9995 against -0.4362). Node2vec is the only embedding of the three
whose ranks vary by site to any meaningful extent, 8 to 10 rank positions
against spectral's 1. Whichever embedding the project ends up using, the
reported reason for preferring spectral does not survive.

---

## 6. What this does not explain, and what is still hypothesis

Stated separately so it is not read as settled.

- ~~**Whether fixing the embedding changes the predictions is untested.**~~
  **TESTED 2026-10-27, and the answer is no.** Four runs were authorised and
  completed: LOO for Generic and Liver FC, predict-all for Control and Liver
  FC. The corrected embedding has 0 of 32 degree dimensions and produces the
  same median rank, a slightly worse paired result on Generic, and worse
  site-invariance. Section 9.
- **Option D is probably the wrong lever, on this evidence.** Dropping the K and
  S blocks cannot recover site specificity that is not in the embedding to
  begin with; Section 4 says it is not. The K-block mechanism in `CATCH_UP.md`
  (c) is not disproved as an *additional* contributor, but it is not the
  binding constraint.
- **Whether the differential-site results are affected: no.** They use the
  measurements directly and never touch the graph, the embedding or the model.
- **Whether a better embedding would predict better is also untested.** Lower
  site rho means sites are distinguishable, which is necessary for
  site-specific prediction and not sufficient. Nothing here shows that any
  filter in Section 4.1 improves held-out ranks.
- **The pathway cliques are implicated but not proven causal.** Section 3.1
  shows kinase degree is 95.5% `protein_same_pathway`, so the ranking is
  pathway-annotation breadth. Whether capping or reweighting those cliques
  restores the spectral gap was not tested; only correlation-edge filters were.
- **Node2vec's apparent site-awareness is not established.** The rank SD
  comparison in Section 2 is single-seed on both sides, and findings 11 and 12
  put node2vec's seed-to-seed variation above the effect sizes this project
  reports. Its 8 to 10 position spread could be walk resampling rather than
  site information. The spectral side of the comparison is unaffected, since
  that embedding is deterministic given the graph.
- **The kinase-side degeneracy is unexplained.** Several filters drive kinase
  rho to 0.92-0.94 while fixing the site side. Why pruning correlation edges
  between *sites* collapses the *kinase* cloud is not established here.

---

## 7. Remedies, cheapest first

1. ~~**Discard the whole null space, not one eigenvector.**~~ **DONE, and it
   does not work.** Implemented on branch `fix/spectral-nullspace` as an opt-in
   `drop_null_space` flag, with the default left at the shipped behaviour so
   existing results stay reproducible. It removes the degree *dimensions*
   (17 of 32 down to 0 of 32) and not the degree *channel* (rho -0.9995 to
   -0.9797), and the median held-out rank is unchanged. Section 9. Keep the
   branch, because the bug is real and the patch is correct, but do not expect
   it to buy anything on its own.

   **1a. Row-normalize the embedding instead.** This is the successor
   candidate, the same size of change, and it targets what Section 9 measures:
   scale each node's embedding vector to unit length, the Ng-Jordan-Weiss
   convention, which removes the sqrt(degree) norm by construction. UNTESTED on
   ranks. About 1.6 h per LOO run at the measured rate.
2. **Report degree as a baseline.** Ranking kinases by degree alone is now
   known to reproduce the model at rho = -0.9995. It belongs alongside the
   popularity baseline in any results table, because it bounds what the model
   adds over connectivity, which on this evidence is close to nothing.
3. **Tighten the correlation-edge magnitude threshold.** On Section 4.1, move
   from the pooled 80th percentile to roughly the 95th, which is the top 25% of
   the current pool. It gives the largest spectral gap of any variant that
   keeps all 3,228 sites, takes site rho from 0.8882 to 0.3685, has the best
   kinase rho of any variant, and adds only 5 components. Prefer it to a
   sample-count floor, which fragments the graph, and do not use FDR control,
   which changes almost nothing. This is a change to the graph build in `src/`.
4. **Reconsider the embedding choice.** Phase 2's preference for spectral over
   node2vec rests on a comparison that the degree channel explains. Once 1 is
   done the comparison should be redone, since node2vec is currently the only
   embedding producing site-varying ranks. Redo it over several seeds, because
   the rank SD evidence for node2vec is single-seed and may be walk noise.
5. **Then, and only then, rerun.** Re-running the pipeline before 1 and 3 would
   re-measure the degree channel.

---

---

## 8. Reproducing this

All read-only, writing only into `thesis/analysis/out/relational_only/`.

```
.venv312/Scripts/python -m thesis.analysis.embedding_geometry
.venv312/Scripts/python -m thesis.analysis.spectral_spectrum
.venv312/Scripts/python -m thesis.analysis.spectral_nullspace
.venv312/Scripts/python -m thesis.analysis.spectral_fix --k 72
.venv312/Scripts/python -m thesis.analysis.degree_rank
.venv312/Scripts/python -m thesis.analysis.corr_edge_variants --k 96
```

The Section 9 runs are NOT read-only in the same sense: they execute the
models on branch `fix/spectral-nullspace`, though they still write only into
`thesis/analysis/out/spectral_fixed/`. About 4.8 h in total at the measured
rate.

```
git checkout fix/spectral-nullspace
.venv312/Scripts/python -m thesis.analysis.test_spectral_nullspace_fix
.venv312/Scripts/python -m thesis.analysis.spectral_fixed_runs --mode estimate
.venv312/Scripts/python -m thesis.analysis.spectral_fixed_runs --mode infer --network control
.venv312/Scripts/python -m thesis.analysis.spectral_fixed_runs --mode infer --network liver_fc
.venv312/Scripts/python -m thesis.analysis.spectral_fixed_runs --mode loo   --network generic
.venv312/Scripts/python -m thesis.analysis.spectral_fixed_runs --mode loo   --network liver_fc
.venv312/Scripts/python -m thesis.analysis.spectral_fixed_compare
.venv312/Scripts/python -m thesis.analysis.spectral_fixed_whynot
```

Recorded as findings 27 to 37 and caveat 12 in `STATE.md`.

---

## 9. The fix was run, and the diagnosis does not survive it

Added 2026-10-27, after four authorised runs on branch `fix/spectral-nullspace`
with `drop_null_space=True` and everything else identical: the frozen 6,909
trials, 415 scored candidates, the original negative sampler, spectral d=32,
`random_state=42`. Every leave-one-out contrast below is paired on
`trial_index` over the 6,894 trials rankable in all runs.

### 9.1 The embedding changed exactly as predicted

Degree dimensions, as R^2 > 0.99 against sqrt(degree) over giant-component
nodes:

| Network | Shipped | Fixed |
|---|---|---|
| Control | 17 / 32 | **0 / 32** |
| Liver FC | 16 / 32 | **0 / 32** |

### 9.2 The ranking did not

| Quantity | Control shipped | Control fixed | Liver FC shipped | Liver FC fixed |
|---|---|---|---|---|
| degree vs mean rank, rho | -0.9995 | **-0.9797** | -0.9973 | **-0.9802** |
| eta^2, rank by kinase | 0.9997 | **0.9999** | 0.9997 | **1.0000** |
| median per-kinase rank SD | 1.05 | **0.52** | 1.18 | **0.24** |
| kinases with SD <= 1 rank | 48.2% | **76.6%** | 42.2% | **93.7%** |

Held-out ranks, leave-one-out:

| Network | Median, shipped | Median, fixed | Median paired diff | Rank-biserial | Favours |
|---|---|---|---|---|---|
| Generic | 79.0 | 79.0 | -2.0 | -0.625 | shipped |
| Liver FC | 74.0 | 74.0 | 0.0 | -0.273 | tie, 4,499 / 6,894 identical |

So the fix removed the dimensions that *are* degree, left the ranking
degree-ordered, did not improve the median rank, came out slightly worse on
Generic, and made the ranking **more** site-invariant rather than less.

### 9.3 Why, measured

The degree signal is not confined to the null space. For the normalized
Laplacian, `L_norm u = lambda u` is equivalent to
`L_rw (D^-1/2 u) = lambda (D^-1/2 u)` where `L_rw = D^-1 (D - A)`, so every
eigenvector of `L_norm` is `u = D^1/2 v` for an eigenvector `v` of the
random-walk Laplacian. Entry `i` of **any** such eigenvector carries a factor
`sqrt(deg_i)`. The null space is only the extreme case, where `v` is constant
on a component and `u` is exactly sqrt(degree).

What survives the fix is therefore the per-node **scale**: node `i`'s embedding
has norm `sqrt(deg_i) * ||v_i||`. Measured Spearman between a kinase's
embedding-vector norm and its degree:

| Network | Shipped | Fixed |
|---|---|---|
| Control | +0.9795 | **+0.9438** |
| Liver FC | +0.9444 | **+0.9444** |

Pair features `[K, S, |K-S|, K*S]` are built from those vectors, so a
high-degree kinase still enters the classifier with a larger vector and keeps
its global position. Scaling every row to unit length removes the norm by
construction, which makes it the next candidate. **Whether that changes
held-out ranks is untested.**

### 9.4 What this costs the rest of the document

- Sections 3 and 3.1 stand as statements about the embedding. The causal
  reading of them, null space to degree to rank, is **refuted** as the
  explanation of the ranking: the arrow from the null space to the rank is not
  there. The arrow from degree to rank remains, at rho = -0.98.
- Section 4's site collapse is untouched by this. It was never attributed to
  the null space, and remedy 3, the correlation-edge threshold, is unaffected
  and is now the more promising of the two graph-side changes.
- Section 5's account of findings 9, 13, 14 and 26 still works through
  *degree*, which is still the ranking key. It no longer works through the
  *bug*.
- Finding 35's direction is worth noting against the random-negatives result:
  two independent corrections, each well motivated and each targeting a
  different real defect, have both **increased** site-invariance (eta^2 to
  1.0000 in both cases). That makes a shared downstream cause more likely than
  a coincidence, and the pair features are the obvious place to look next.

### 9.5 Honest note on the measurement

My first version of the Section 9.1 table measured R^2 over all nodes instead
of giant-component nodes and reported 0 of 32 for the shipped setting too,
which would have made the bug look nonexistent. ARPACK returns an arbitrary
basis of the null space, so a returned vector is a mixture of the
per-component degree vectors with a different scale on each; over all nodes at
once such a mixture does not fit sqrt(degree) even though each
component-restricted piece fits it exactly. The giant-component scope is the
one `spectral_fix.py` used and the one reported here.
