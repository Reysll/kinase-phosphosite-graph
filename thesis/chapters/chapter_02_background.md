# CHAPTER II: BACKGROUND AND RELATED WORK

> Draft. Citation verification completed 2026-10-10. `[VERIFIED]` = DOI or PMC
> link confirmed, recorded in `thesis/references.md`. `[PARTIAL]` = freely
> accessible but no DOI exists. Two items remain unresolved and are marked
> inline: the protein-interaction source file has no recorded provenance, and
> the liver dataset's attribution to Gao et al. is blocked by a cohort-size
> mismatch. No citation is guessed.

## 2.1 Protein Phosphorylation and Kinase Signaling

Protein phosphorylation is the reversible attachment of a phosphate group to a
specific amino acid residue, most often serine, threonine, or tyrosine. The
reaction is catalyzed by protein kinases and reversed by protein phosphatases.
Because adding a phosphate changes the local charge and shape of a protein, it
can switch enzymatic activity on or off, create or destroy a binding surface for
another protein, or mark a protein for degradation. Phosphorylation is therefore
one of the principal mechanisms by which cells transmit information from the
membrane to the nucleus and by which they respond to stress, growth signals, and
damage.

Two properties of this system matter for the computational problem addressed in
this thesis. First, phosphorylation is **site-specific**: a single protein may
carry many phosphorylation sites, and different sites on the same protein can be
controlled by different kinases and produce different functional outcomes. The
prediction target is therefore a kinase-to-site relationship, not a
kinase-to-protein relationship. Second, phosphorylation is **context-dependent
and transient**: which kinase acts on a given site depends on the tissue, the
cell state, and the presence of disease. A kinase-substrate relationship
observed in one cellular context does not necessarily hold in another.

The set of human protein kinases, collectively the kinome, comprises roughly 500
enzymes. Dysregulation of kinase signaling is a recurring feature of cancer, and
kinases are among the most successfully drugged protein families, which makes
the identification of disease-relevant kinase-substrate relationships
translationally important as well as biologically interesting.

## 2.2 The Annotation Gap and the Streetlight Effect

Although high-throughput mass spectrometry has catalogued hundreds of thousands
of human phosphorylation sites, the enzyme responsible for most of them is
unknown. Deznabi et al. (2020) state that "more than 95% of the reported human
phosphosites have no known kinases" [VERIFIED]. The gap is not merely large, it
is also **structurally uneven**. Maier et al. (2025) report that approximately
150 kinases have no known substrates at all, and that roughly 80% of annotated
phosphosites are attributed to just 20% of kinases [VERIFIED].

This imbalance produces a self-reinforcing bias often called the streetlight
effect: well-studied kinases accumulate more annotations, which makes them
easier to study further, while understudied kinases remain dark. For machine
learning the consequence is direct. A model trained on existing annotations
inherits the imbalance, and an evaluation built from the same annotations
rewards a model for performing well on kinases that were already well
characterized. Chapter V returns to this point with a quantitative
demonstration: in the evaluation set used here, 85.2% of trials involve
well-annotated kinases, and a baseline that ranks candidates purely by how many
substrates they already have is a surprisingly strong competitor.

## 2.3 Resources for Kinase-Substrate Annotation

Several curated resources supply the evidence integrated in this work.

**PhosphoSitePlus** is the principal source of experimentally determined
kinase-substrate associations and of phosphorylation site definitions
(Hornbeck et al., 2012) [VERIFIED]. Its curated kinase-substrate table provides
the positive labels used throughout this thesis.

**DEPOD** catalogues human phosphatase-substrate relationships, supplying the
dephosphorylation counterpart to the kinase annotations (Duan et al., 2015)
[VERIFIED].

**MSigDB** and **PTMsigDB** provide pathway membership at the protein level and
at the post-translational-modification site level respectively, allowing
proteins or sites that participate in the same biological process to be linked
(Liberzon et al., 2015; Krug et al., 2019) [VERIFIED].

**Protein-protein interaction data** supplies physical and high-confidence
functional associations between proteins. The specific source of
`high_confidence_score.csv` is not recorded anywhere in the repository and is
**not yet attributable**; see `thesis/references.md`.

**PTMcode2** provides evidence of co-evolution between modification sites, which
captures functional coupling not visible in interaction or pathway data
(Minguez et al., 2015) [VERIFIED].

These resources are complementary rather than redundant: each supplies a
different kind of evidence about how proteins and sites relate, which is the
premise of the heterogeneous graph described in Chapter III.

## 2.4 Computational Approaches to Kinase-Substrate Prediction

Existing methods for predicting which kinase acts on which site fall broadly
into three families.

### 2.4.1 Sequence-Based Methods

The earliest and still most widely used approaches learn the amino acid motif
surrounding sites phosphorylated by a given kinase, then scan candidate sites
for matches to that motif. Motif-based scoring is fast, interpretable, and
effective for kinases with sharply defined specificity. Its limitations are
twofold. Many kinases share highly similar motifs, so sequence alone cannot
distinguish family members, and a motif model requires enough known substrates
to be estimated in the first place, which excludes exactly the understudied
kinases that matter most.

### 2.4.2 Structure-Based Methods

Structure-aware approaches model the physical compatibility between a kinase
catalytic domain and a candidate substrate peptide. These methods can in
principle generalize to kinases with few known substrates, since they rely on
structural rather than statistical evidence. They are limited by the
availability and accuracy of structural models and by the computational cost of
modelling every kinase-site pair.

### 2.4.3 Network-Based Methods

Network methods represent proteins, sites, and their relationships as a graph
and predict kinase-substrate links from graph structure. The rationale is that
a kinase is more likely to act on a site if the two are close in a network of
pathway membership, physical interaction, co-expression, or co-evolution, even
when sequence similarity is uninformative. This family is the one extended in
this thesis. Its appeal is that it naturally accommodates heterogeneous evidence
and can be made context-specific by changing the graph. Its risk, documented in
Chapter V, is that graph structure encodes connectivity as well as biology, and
a model may learn the former while appearing to learn the latter.

Methods that explicitly target understudied kinases have also been proposed.
DeepKinZero (Deznabi et al., 2020) uses zero-shot learning to transfer
information to kinases with no training substrates [VERIFIED], an approach
motivated by the same annotation imbalance that motivates the stratified
modeling in this thesis.

## 2.5 Graph Representation Learning

Converting a graph into fixed-length vectors suitable for a downstream
classifier is the central modelling step in network-based prediction. Two
families are compared in this work.

### 2.5.1 Random-Walk Embeddings

DeepWalk (Perozzi et al., 2014) and node2vec (Grover & Leskovec, 2016)
[both VERIFIED]
generate sequences of nodes by random walks over the graph and then apply a word
embedding model, treating node sequences as sentences. Nodes that frequently
co-occur in walks receive similar vectors. node2vec adds two parameters, `p` and
`q`, that bias walks toward breadth-first or depth-first exploration; with
`p = q = 1` the procedure reduces to DeepWalk. Random-walk embeddings capture
local neighbourhood structure well and scale to large graphs.

Their central weakness for this application is **degree sensitivity**.
High-degree nodes appear in a disproportionate share of random walks, so their
embeddings become central in the resulting vector space. In a kinase-substrate
graph the highest-degree kinases are precisely the best-annotated ones, so the
embedding amplifies the annotation bias the method is intended to overcome.
Chapter V documents this effect empirically and also shows that random-walk
embeddings carry a second, less widely discussed problem: substantial
instability with respect to the random seed.

### 2.5.2 Spectral Embeddings

Spectral embedding derives coordinates from the eigenvectors of a graph
Laplacian. Using the symmetric normalized Laplacian

$$L = I - D^{-1/2} A D^{-1/2}$$

where `A` is the adjacency matrix and `D` the diagonal degree matrix, each row
and column of `A` is scaled by the inverse square root of the corresponding
node's degree. High-degree nodes are therefore down-weighted by construction,
which makes spectral embedding a natural remedy for the degree sensitivity of
random-walk methods. The eigenvectors associated with the smallest non-trivial
eigenvalues encode the coarse community structure of the graph. Spectral methods
have a long history in clustering and dimensionality reduction (von Luxburg, 2007
[VERIFIED]; Shi & Malik, 2000 [VERIFIED]; Ng, Jordan, & Weiss, 2002 [PARTIAL:
no DOI, NeurIPS-era proceedings]).

Two properties relevant to later chapters follow from this definition. Spectral
embedding is **deterministic** given the graph and the number of components
requested, so it has no seed-stability problem. And because only the leading
eigenvectors are retained, it represents **global** community structure rather
than local neighbourhood detail, which turns out to determine what kind of
context signal the method can and cannot see.

### 2.5.3 Degree Bias as a Known Failure Mode

Degree bias in graph embeddings is a recognized issue, and the standard
diagnostic is to compare a model against a null graph in which degrees are
preserved but edges are rewired. If performance survives randomization, the
model is responding to the degree sequence rather than to the specific pattern
of connections. This thesis applies that diagnostic directly in Chapter V, using
double-edge-swap randomization that preserves each node's degree exactly.

## 2.6 Context-Specific Biological Networks

A generic interaction network represents an average over cell types and
conditions. Because phosphorylation is context-dependent, a generic network may
mislead when applied to a particular tissue or disease. Context-specific network
construction addresses this by reweighting or filtering a generic network using
condition-specific measurements, most often transcriptomic or proteomic
abundance data.

A common approach is to add edges between molecular entities whose abundances
co-vary across samples from the condition of interest, on the assumption that
co-variation reflects shared regulation. This is the strategy adopted in
Chapter III, where phosphosite-level correlations computed from liver cancer
phosphoproteomics are added to the generic graph as new edges. The approach has
a known vulnerability that Chapter VI discusses at length: correlation edges add
connectivity as well as information, and if a model is sensitive to connectivity
alone then adding edges can appear to help for reasons unrelated to the biology
being encoded.

The liver cancer measurements used here come from an integrated proteogenomic
characterization of hepatitis-B-related hepatocellular carcinoma (Gao et al., 2019)
[VERIFIED publication, PROVENANCE UNRESOLVED], which provides matched protein
and phosphoprotein abundance for tumor and adjacent non-tumor tissue. Note that
Gao et al. report 159 patients whereas the file used here contains 18 tumor and
18 control samples; the attribution is therefore provisional pending
confirmation (see `thesis/references.md`).

## 2.7 Inferring Kinase Activity from Phosphoproteomics

A complementary question to kinase-substrate prediction is kinase **activity**
inference: given measured changes in phosphorylation across conditions, which
kinases appear more or less active? Kinase-Substrate Enrichment Analysis (KSEA)
answers this by testing whether the known substrates of a kinase are
collectively shifted relative to the global distribution of changes (Wiredja
et al., 2017) [VERIFIED]. The score takes the form

$$\text{score} = \frac{(\bar{s} - \bar{p})\sqrt{m}}{\delta}$$

where `s̄` is the mean log fold-change of the kinase's substrates, `p̄` the
global mean, `m` the number of substrates observed, and `δ` the global standard
deviation. The `√m` term in the numerator rewards larger substrate sets through
a smaller standard error.

KSEA and network-based ranking answer different questions. Ranking asks which
kinase is structurally most plausible for a given site; KSEA asks which kinases
show coordinated substrate-level change. Chapter V reports both and finds that
they nominate largely disjoint kinase sets, which is informative about what each
method measures.

## 2.8 Evaluation Practice and Its Pitfalls

Three aspects of evaluation design recur in this thesis.

**Leave-one-out evaluation.** Holding out a single known kinase-substrate edge
and asking the model to recover it is the standard protocol. It uses the
available annotations efficiently, which matters when positives are scarce.

**Leakage through representation.** When a graph embedding is computed on a
graph that still contains the edges being predicted, information about the
answer propagates into the node vectors even if the specific edge is removed
from the classifier's training set. Preventing this requires removing the target
edge type from the graph before embedding, not merely from the label set. The
implementation in this thesis does so, and Chapter IV documents the effect.

**Baselines.** A model that outperforms random guessing has not necessarily
learned anything useful, because the relevant comparison is against the simplest
method that exploits the same statistical structure. In a domain with severe
annotation imbalance, the natural trivial baseline is to rank candidates by
annotation count. Chapter V shows that this baseline is strong, and that
evaluating against it changes the interpretation of the headline results
substantially.

## 2.9 Summary and the Gap Addressed

Existing kinase-substrate prediction methods are predominantly context-free and
are biased toward well-annotated kinases. Network-based methods offer a route to
incorporating disease context by modifying graph structure, and spectral
embeddings offer a principled way to reduce the degree sensitivity of
random-walk methods. What has been less thoroughly examined is whether the
apparent benefit of adding context to a biological graph reflects the biological
content of the added edges or merely their effect on connectivity, and whether
such models outperform the annotation prior they inherit.

This thesis builds a heterogeneous protein-phosphosite network, extends it with
liver cancer phosphoproteomic context, evaluates kinase-substrate prediction
across four network variants and multiple embedding strategies, and then
subjects the resulting claims to null-model, seed-stability, site-dependence,
and baseline testing. Chapter III describes the data and graph construction.
