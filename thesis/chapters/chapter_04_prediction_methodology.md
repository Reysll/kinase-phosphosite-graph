# CHAPTER IV: PREDICTION METHODOLOGY

This chapter describes how the graphs of Chapter III are turned into kinase
rankings, and how those rankings are evaluated. It covers methods only. No
result is reported here; Chapter V presents results.

The description separates two things that are easy to conflate. Sections 4.2
through 4.5 describe the predictor, meaning how a numerical score is assigned to
a kinase and phosphosite pair. Sections 4.6 through 4.9 describe the evaluation
apparatus, meaning how those scores are turned into a measurement and what
controls establish whether the measurement means anything.

## 4.1 Overview

The prediction task is framed as scoring pairs. Given a phosphosite and a set of
candidate kinases, the method assigns each candidate a score and ranks them.
A good method ranks the kinase that is actually responsible near the top.

The pipeline has four stages:

1. **Embed.** Every node in the graph is mapped to a fixed-length vector that
   summarizes its position in the network.
2. **Combine.** A kinase vector and a site vector are combined into a single
   feature vector describing the pair.
3. **Classify.** A logistic regression model trained on known kinase-substrate
   edges assigns each pair a probability.
4. **Rank.** Candidates are sorted by score, and the position of the correct
   kinase is recorded.

The design keeps these stages independent. Embeddings are computed once per run
and reused across all trials, the pair construction is a fixed arithmetic
operation, and the classifier is retrained per trial. This separation is what
makes it possible to swap one embedding method for another and attribute any
change in performance to that substitution alone.

## 4.2 Node Embeddings

### 4.2.1 The Strategy Interface

Embedding methods are interchangeable through an abstract interface. A strategy
implements a single operation, mapping a graph to a dictionary of node
identifiers and vectors:

$$\text{fit}: G \rightarrow \{v \mapsto \mathbf{x}_v \in \mathbb{R}^d\}$$

Each strategy also declares whether it requires a directed or an undirected view
of the graph, which the runner honors when constructing the embedding graph.
Adding a new method requires implementing this one operation and nothing else.
Three strategies are used in this work.

### 4.2.2 Random-Walk Embeddings

The first strategy is node2vec, which learns vectors by sampling random walks
from each node and training a skip-gram model so that nodes appearing in similar
walk contexts receive similar vectors. It operates on the directed graph.
Because the return and in-out parameters are both left at one, the walks are
unbiased and the method reduces to DeepWalk, a deliberate choice that avoids
introducing two further hyperparameters whose tuning could not be justified on
this data.

Random-walk methods have a known structural property that matters here. Walks
visit high-degree nodes disproportionately often, so hub nodes appear in many
contexts and acquire vectors that are, in effect, more central in the embedding
space. Section 2.5.3 discusses this bias in the literature. It is relevant
because the graph described in Chapter III concentrates 90.6% of its edges in
pathway co-membership, which is densest around well-studied proteins.

### 4.2.3 Spectral Embeddings

The second strategy embeds nodes using the eigenvectors of the symmetric
normalized Laplacian,

$$L = I - D^{-1/2} A D^{-1/2}$$

where $A$ is the adjacency matrix and $D$ is the diagonal degree matrix. The
$D^{-1/2}$ terms on both sides are what distinguish this from the unnormalized
Laplacian, and they are the reason the method is less dominated by degree: each
node's row and column are scaled by its own connectivity, so a hub does not
simply outweigh everything else by virtue of having more edges.

The implementation requests the smallest eigenvalues by magnitude using a sparse
iterative eigensolver, discards the trivial eigenvector associated with the
near-zero eigenvalue, and retains the next $d$ eigenvectors as coordinates. The
normalized Laplacian is defined only for symmetric adjacency, so this strategy
declares that it requires an undirected graph and directed edges are symmetrized
before the decomposition.

One practical constraint follows from the eigensolver. If the graph fragments
into many disconnected components, the near-zero eigenvalue acquires high
multiplicity and the decomposition becomes unstable. This is why the generic
network is built from the dedicated generic graph file rather than from the
liver graph restricted to base relations, a substitution that would isolate
every liver site whose only connections are the correlation edges being removed.

### 4.2.4 Concatenation

The third strategy concatenates the two, producing vectors of dimension $2d$.
It exists because the two methods address different aspects of the degree
problem, spectral reducing the global dominance of hubs while leaving
within-group ordering largely intact. Each sub-strategy receives a graph of the
directedness it requires rather than sharing one view.

### 4.2.5 Dimensions and Settings

Node2vec is run with walk length 10, 25 walks per node and a context window of
5. Local leave-one-out runs and the predict-all pass use $d = 32$ per
sub-embedding, giving 64-dimensional nodes under concatenation. Earlier cluster
runs used $d = 64$, so feature counts quoted from those runs are not comparable
and are not used in this thesis. Node2vec is run with four worker threads and a
batch size of 10,000 words in all production runs. These differ from the class
defaults, so any new caller must pass them explicitly to remain comparable.

## 4.3 The Leakage Guard

The most consequential design decision in the pipeline concerns what the
embeddings are allowed to see.

The prediction target is the `phosphorylates` relation, which encodes known
kinase-substrate associations. If those edges were present in the graph when the
embeddings were computed, a kinase and its substrate site would be adjacent, the
embedding would place them close together, and the classifier would recover that
proximity. The resulting score would measure whether an edge had been memorized,
not whether it could be predicted.

The embedding graph is therefore constructed with **all `phosphorylates` edges
removed**. Embeddings see protein-site membership, dephosphorylation, pathway
co-membership, protein interaction, co-evolution and, where applicable,
correlation edges, but never a kinase-substrate association. Those edges are
reintroduced only as training labels for the classifier.

A direct consequence is that embeddings depend on the network and the strategy
but not on the trial. They are computed once per run, before the trial loop
begins, and reused across every trial. Moving the embedding computation inside
the loop would be both wasteful and, if the held-out edge were restored to the
graph first, a reintroduction of the leak.

## 4.4 Pair Features

A kinase vector $\mathbf{k}$ and a site vector $\mathbf{s}$, each of dimension
$d$, are combined into a single feature vector of dimension $4d$:

$$\mathbf{f}(\mathbf{k}, \mathbf{s}) = [\, \mathbf{k} \;\|\; \mathbf{s} \;\|\; |\mathbf{k} - \mathbf{s}| \;\|\; \mathbf{k} \odot \mathbf{s} \,]$$

where $\|$ denotes concatenation and $\odot$ the elementwise product. The four
blocks carry different information. The first two preserve the identity of each
partner, allowing the model to learn that some kinases are more often correct
and some sites more often targeted. The absolute difference measures
dissimilarity per dimension and is symmetric. The elementwise product measures
agreement per dimension and is the term that can express interaction between the
two vectors rather than properties of either alone.

At $d = 32$ this yields 128 features. Pairs in which either partner lacks an
embedding are dropped rather than imputed, which is the mechanism by which five
of the 420 candidate kinases are never scored.

The first two blocks deserve a note, because they are also a vulnerability.
Since $\mathbf{k}$ enters the feature vector directly, the model can achieve
nontrivial training accuracy using kinase identity alone, without reference to
the site. Section 4.9.3 describes the diagnostic designed to detect exactly
this.

## 4.5 The Classifier

Scoring uses $L_2$-regularized logistic regression with the liblinear solver, a
1,000-iteration cap and balanced class weights.

Positives are the known kinase-substrate edges in the training set. Negatives
are drawn per site from the candidate kinases not known to act on that site,
capped at 50 per site, and balanced class weighting compensates for the
resulting imbalance by reweighting the classes in inverse proportion to their
frequency.

The selection of those 50 requires a precise statement, because the
implementation does not do what the word "sample" would suggest. The code
filters the candidate list to the kinases not known to act on the site and then
takes the **first 50 entries of the resulting list**. No randomization is
applied. The candidate list is itself constructed with a sort on the node
identifier, so its order is **alphabetical by gene symbol**.

The consequence is that negatives are not a sample of the candidate space but a
fixed alphabetical prefix of it. Measured over all 2,459 training sites, only
**55 of the 420 candidate kinases are ever used as a negative example**, and
each of those appears at essentially every site. **361 of the 415 scored
kinases, 87.0%, never appear as a negative at any site.** The decision boundary
is therefore fitted against an alphabetically determined subset of the kinases
the model is later asked to rank.

This is a defect in the training design and is reported as a limitation.
Whether it biased the results was tested directly, by comparing the rank
distributions of never-negative and ever-negative kinases in the existing
output; Section 5.9 gives the figures. No comparison reached significance, so
the defect is documented without any accompanying claim that it distorted the
rankings. One consideration also affects how much room there is for such an
effect: as Section 4.9.3 describes, kinase identity accounts for almost all the
variance in predicted rank, which leaves little that the composition of the
negative set could move.

A second consequence is that negative selection is **deterministic**. Given the
same training positives it returns the same pairs every time, so it contributes
nothing to run-to-run variability.

Treating unobserved pairs as negatives is an assumption, and given the
annotation gap documented in Section 2.2 it is certainly wrong for some pairs.
A kinase with no recorded association to a site may genuinely not act on it, or
may simply never have been tested. The evaluation inherits this limitation, and
Section 4.9.4 describes the control that measures how much of the apparent
performance is attributable to annotation frequency rather than biology.

The choice of a linear model over a more expressive one is deliberate. The
scientific question is whether adding context edges changes the information
available in the graph. A linear classifier makes the comparison interpretable
by holding model capacity fixed, so that a change in measured performance is
attributable to the representation rather than to the classifier's ability to
exploit it.

## 4.6 Leave-One-Out Evaluation

### 4.6.1 The Frozen Trial Set

Evaluation holds out one known kinase-substrate edge at a time, retrains, and
measures where the held-out kinase lands in the ranking for its site.

Trials are restricted to **multi-kinase sites**, meaning sites with more than one
annotated kinase. On a single-kinase site, removing the only known association
leaves no positive label for that site and the adjusted rank of Section 4.6.3
becomes undefined.

The trial set is **frozen**, meaning the held-out edges are selected once, saved,
and reused identically by every run. This is what makes the four networks and
the three embedding strategies directly comparable: each is evaluated on exactly
the same trials in exactly the same order, so differences cannot arise from
sampling. The frozen set contains **6,909 trials across 2,459 sites**.

Of those, **15 trials are unrankable in every run**, because the held-out kinase
lacks an embedding and therefore receives no score. These are excluded
consistently rather than imputed, so all metrics in this thesis are computed on
**n = 6,894**. Of the 420 defined candidate kinases, **415 are scored** in every
run, the remaining five lacking embeddings for the same reason.

### 4.6.2 A Single Trial

Each trial proceeds identically:

1. Remove the held-out edge from the training positives.
2. Train the classifier on the remaining positives and on the negatives selected
   as described in Section 4.5. Those negatives are deterministic, so they are
   identical across trials except at the site whose true-kinase set the holdout
   has just changed.
3. Score all 415 candidate kinases against the held-out site.
4. Sort by score, highest first.
5. Record the position of the held-out kinase.

Embeddings are not recomputed, for the reason given in Section 4.3.

### 4.6.3 Adjusted Rank

The naive rank of the held-out kinase is misleading on multi-kinase sites.
Suppose a site has known kinases $k_1$ and $k_2$, and $k_1$ is held out. The
edge for $k_2$ remains in training, so the model has been told directly that
$k_2$ acts on this site and will rank it highly. If $k_2$ ranks above $k_1$, it
pushes $k_1$ down by one position for a reason that has nothing to do with the
quality of the prediction for $k_1$.

The **adjusted held-out rank** removes this artifact. All other known true
kinases for the site are dropped from the scored list before the held-out
kinase's position is measured. If the scored order is $[k_7, k_2, k_4, k_1]$
then the naive rank of $k_1$ is four, and dropping $k_2$ gives an adjusted rank
of three.

This is the primary metric throughout this thesis. The unadjusted rank is
retained in the output for backward compatibility but is not used for any
comparison.

### 4.6.4 Within-Category Rank

Candidate kinases are assigned to one of three categories by their number of
annotated substrate sites in PhosphoSitePlus, counted over the full generic
graph rather than the evaluation subset:

| Category | Substrate sites | Kinases |
|---|---|---|
| poor | at most 5 | 166 |
| average | 6 to 20 | 133 |
| rich | more than 20 | 121 |

The counts above are over the 420 defined candidates. Restricted to the 415
actually scored, the split is 162 poor, 132 average and 121 rich.

The motivation is a fairness problem in the global ranking. A kinase with three
annotated substrates and one with three hundred do not compete on equal terms,
because the model has seen a hundred times more positive examples for the
second. A poorly characterized kinase can therefore rank low for reasons of data
volume rather than biology, which is a serious objection given that
understudied kinases are precisely the ones this work aims to say something
about.

The **within-category adjusted rank** addresses this by ranking each kinase only
against others in its own category. It mirrors the adjusted global rank exactly:
other true kinases in the same category as the held-out kinase are dropped
before its position is measured. The denominator falls from 415 to roughly 120
to 166 depending on the category, so within-category ranks are not comparable to
global ranks in absolute terms and the two are always reported separately.

## 4.7 Metrics

Three summaries are reported per run.

**Median adjusted rank.** The median over trials of the adjusted held-out rank.
The median is preferred to the mean because the rank distribution is skewed and
bounded below at one.

**Hit@10 and Hit@20.** The fraction of trials in which the adjusted rank is at
most 10 or at most 20 respectively. These are absolute cutoffs, not percentiles.
Earlier project documentation labelled them "Top-10%" and "Top-20%", which
wrongly implies a proportional cutoff; that naming is not used here.

Because all runs share the frozen trial set, every metric is paired at the trial
level, which is what permits the paired tests described in Chapter V.

## 4.8 Predict-All Inference

Leave-one-out answers a narrow question: for an edge known to exist, how highly
is the correct kinase ranked. It cannot produce a rank for every kinase at every
site, because it only ever evaluates the kinase that was held out, and only on
the multi-kinase sites in the frozen set. Comparing how a given kinase's ranking
changes between networks requires the full matrix.

The predict-all pass produces it, scoring every candidate kinase against every
one of the 3,228 fold-change-valid liver sites, for each network and embedding
combination, with no leave-one-out holdout.

Avoiding leakage in this setting requires **two-tier scoring**, because most
pairs were never positive labels and a few were:

- Pairs that are **known true edges**, approximately 736 within the 3,228 sites,
  are scored by a model trained with **only that single edge masked**. Each such
  pair gets its own model.
- **All other pairs** are scored by one global model trained on all known edges.
  Nothing can leak for a pair that was never a positive label.

This is far cheaper than it appears, since the expensive masked refits number in
the hundreds rather than the thousands.

Output is in long form, one row per site and kinase pair, carrying the score,
the rank within that site and a flag for whether the pair is a known true edge.

One asymmetry must be stated. The generic network covers only **391 of the 3,228
sites**, about 12%, because its site population is essentially the
PhosphoSitePlus-annotated set rather than the broader liver-detected set. Its
site population is therefore a different population, and the generic network is
reported for reference but never placed inside a contrast with the others.

## 4.9 Evaluation Controls

Measuring that one configuration outranks another establishes only that the
numbers differ. It does not establish that the difference reflects the
biological structure the method is supposed to exploit. Four controls were
implemented to test that interpretation directly. Each targets a specific
alternative explanation. All are read-only with respect to the project's data
and results, and all write to a separate analysis directory.

### 4.9.1 The Degree-Preserving Null

**Alternative explanation: the context edges matter only through how many of
them each node receives, not through which nodes they connect.**

The null replaces the liver correlation edges with random edges whose **degree
sequence is identical**, constructed by repeated double-edge swaps. A swap takes
two edges $(a,b)$ and $(c,d)$ and replaces them with $(a,d)$ and $(c,b)$, which
leaves every node's degree unchanged while destroying the correspondence between
the edge and the measured co-variation that produced it. Swaps are applied at
ten times the edge count, which randomizes essentially every edge.

The randomized network is then run through the identical pipeline. If the real
correlation edges carry information beyond their degree contribution, the real
network should outperform its null. If it does not, then what the method
responds to is connectivity volume, which any random network with the same
degrees would supply.

This control is applied to both the leave-one-out comparison and the predict-all
rank shifts, and separately to within-category ranking.

### 4.9.2 Seed Variance

**Alternative explanation: the difference between networks is smaller than the
run-to-run variation of the method itself.**

This control applies to node2vec and to node2vec only, and the reason is worth
stating because it determines the scope of every seed-related claim in this
thesis.

Node2vec is stochastic in two places, the sampling of random walks and the
initialization of the skip-gram model, so two runs with different seeds produce
different embeddings from the same graph. The spectral strategy is not. It
solves an eigenproblem whose starting vector is drawn from a generator seeded by
an explicit `random_state` parameter, fixed at its default throughout this work,
so for a given graph it returns the same embedding every time. This matches the
determinism noted in Section 2.5.2. Negative selection is deterministic for the
reason given in Section 4.5. Run-to-run variability therefore has exactly one
source, and it is present only under node2vec.

The control accordingly reruns the **node2vec** configuration at several seeds,
holding the network and everything else constant, and measures the spread. The
comparison of interest is between the magnitude of the seed-to-seed variation
and the magnitude of the network-to-network difference. A context effect smaller
than the seed noise is not a context effect.

No seed sweep was run for spectral, because repeating a deterministic
computation would return identical results. The corollary is that seed-based
findings constrain the node2vec results and say nothing about the spectral
ones, which must be qualified by the separate controls instead.

### 4.9.3 Site Dependence

**Alternative explanation: the model applies a fixed ranking of kinases and
barely conditions on the site being scored.**

This is the failure mode Section 4.4 anticipates. If the kinase vector alone
drives the score, a kinase will occupy nearly the same rank at every site, and
the model will be producing a global kinase ordering dressed up as site-specific
prediction.

The diagnostic is a one-way analysis of variance of the predicted rank on kinase
identity, across all sites in the predict-all output, summarized by

$$\eta^2 = \frac{SS_{\text{between}}}{SS_{\text{total}}}$$

which is the proportion of variance in predicted rank attributable to which
kinase it is. A value near one means site identity contributes almost nothing.
The per-kinase standard deviation of rank across sites is reported alongside, as
a directly interpretable companion in units of rank positions.

### 4.9.4 The Popularity Baseline

**Alternative explanation: the apparent performance reflects annotation
frequency, which requires no model at all.**

The baseline ranks candidate kinases by their number of annotated substrate
sites, identically at every site. It uses no graph, no embedding and no
learning.

For the comparison to be fair it must be held to exact parity with the models.
The baseline is therefore constructed from the same frozen edge set the models
train on rather than from the full annotation table, and it scores the same 415
kinases over the same 6,894 rankable trials. An earlier version violated this
parity by drawing counts from the full table, which inadvertently weakened the
baseline; the figures in this thesis use the corrected version.

The baseline is evaluated on the global adjusted rank, within each category, and
stratified by the held-out kinase's category, because a baseline built from
substrate counts is expected to behave very differently for poorly and richly
annotated kinases and a single aggregate figure would obscure that.

## 4.10 Summary

The predictor embeds a graph stripped of its prediction target, combines kinase
and site vectors through concatenation, absolute difference and elementwise
product, and scores pairs with a regularized linear classifier. Evaluation holds
out one edge at a time over a frozen set of 6,909 trials, of which 6,894 are
rankable, ranking 415 candidate kinases and reporting the adjusted rank both
globally and within substrate-count categories. A predict-all pass extends
scoring to every kinase and site pair using two-tier leakage-safe scoring.

Four controls accompany the evaluation: a degree-preserving null, a seed
variance analysis, a site-dependence diagnostic and a popularity baseline. They
were added because each tests a specific way in which a measured difference
could arise without the method working as intended. Chapter V reports what they
found.
