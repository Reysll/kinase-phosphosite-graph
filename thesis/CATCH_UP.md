# Catch-up: where this project stands

Written to be explained out loud. Last updated 2026-10-26.

## (a) What we set out to do

Kinases are enzymes that attach phosphate groups to specific spots on other
proteins. Those spots are called phosphosites, and which kinase hits which site
is the thing we want to predict. The problem is that most of what is known comes
from a general-purpose database that is not specific to any tissue or disease,
and more than 95% of known phosphosites have no kinase assigned to them at all.

Our idea was to make the prediction context-aware. We built a network that holds
proteins and phosphosites as nodes, with different kinds of edge for pathway
membership, protein interaction, co-evolution and so on. Then we added liver
cancer data on top, connecting phosphosites whose measured levels move together
across patients. The hypothesis was that adding liver-specific information would
make the predictions better for liver cancer, and that we could then hand a
short list of kinase and substrate pairs to a wet lab for testing.

## (b) The six things we found

**1. The liver data does not actually add information.** We rebuilt the liver
network with random edges that preserve exactly how many connections each node
has, destroying the real biology but keeping the shape, and the fake network
gives the same answer: 96% of predictions come out with an identical rank. If
randomized biology works as well as real biology, the model is responding to how
connected a node is, not to what it is connected to.

**2. The strongest single result we had was reproduced by that random network.**
We had found a cluster of necroptosis kinases, RIPK1, RIPK3 and MLKL, that
looked more relevant in the cancer networks, and it was the most biologically
exciting thing in the project. The random network reproduces 100% of that shift,
rank for rank, so the cluster is a property of the added connections rather than
of liver cancer.

**3. Changing a random seed moves the result more than changing the biology
does.** When we ran the same configuration with different random seeds, the
seed-to-seed difference was larger than the difference between the generic and
the liver networks, 0.373 against 0.249 on the same effect-size scale. That
means most of what we were reporting as a biological effect was within the noise
of rerunning the same experiment.

**4. The model barely looks at the phosphosite.** We measured how much of the
variation in a kinase's predicted rank comes from which kinase it is, as opposed
to which site it is being scored against, and the answer is 99.97%. The model is
essentially producing one fixed ranking of kinases and repeating it at every
site, which is not site-specific prediction.

**5. A baseline with no model in it beats everything we built.** If you just
rank kinases by how many substrates they already have on record, ignoring the
network entirely, you get a median rank of 25 against 74 for our best model.
There is one real exception and it matters: for poorly studied kinases, which
are the ones we most want to say something about, our models do beat that
baseline by a wide margin.

**6. The differential analysis works, and it is the part we can hand over.**
Comparing phosphosite levels directly between tumour and adjacent healthy tissue
gives 659 sites that change substantially and reliably. This uses only the
measurements, not the network or the model, so none of the problems above apply
to it.

## (c) What the advisor asked for on October 15, and where each stands

**Write it up for a journal, not a conference.** Done as a direction. There is
no deadline, and the robustness checks we ran are exactly what a journal
reviewer would ask for, so they help us rather than hurt us.

**Get the collaborator statistically significant, liver-cancer-associated
substrates.** Done. The 659-site shortlist is in `thesis/advisor_package/`, with
a one-page summary and a CSV.

**Stop using prediction performance as the headline, use comparative kinase
ranking within substrate-count categories instead.** Done, and the answer is
not what we hoped. The category ranking was already built back in June, so we
could test it immediately, and it does not separate the real network from the
random one. The two agree at a correlation of 0.9944.

**Re-examine the preprocessing after the DEPOD count looked wrong.** Done, and
the advisor was right. The file has 115 phosphatases, not 250, and our pipeline
keeps only 65 of those. But that layer is 0.03% of the network, so it cannot
explain any result. The audit found something bigger: 90.6% of the network is
pathway co-membership, which is an annotation rather than a measurement, and
that is very likely why a random network does as well as the real one.

**Explain the correlation sentence that did not make sense.** Done, and written
up in `thesis/GRAPH_WALKTHROUGH.md` for the step-by-step meeting that was
requested. Short version: we do not correlate across all 18 samples. Each pair of
phosphosites is correlated using only the samples where both were detected, and
the most common case is the minimum of 6. Only 8% of our edges use all 18.

## (d) What is still undecided

**What the main result of the paper is.** The ranking work cannot distinguish
itself from a random network, which is a legitimate negative methodological
finding but is not the paper we set out to write. Either the network changes, or
the differential phosphosite work becomes the main contribution and the ranking
work is reported as a cautionary result. This is the decision to make.

**Whether to change the network.** The obvious lever is the pathway edges at
90.6%. Options are capping how large a pathway clique can get, weighting edges
by pathway size, or removing those edges entirely and measuring what is lost.

**Where the liver data came from.** The file holds supplementary tables S4 and
S11 of some hepatocellular carcinoma paper, but we cannot identify which one.
One candidate was ruled out because it has 159 patients and ours has 18 and 18.
Until this is resolved we cannot cite the dataset, and we cannot get the
patient-level clinical information either.

**Whether the samples are matched pairs.** The layout looks like 18 patients
each giving a tumour and a healthy sample, but the file does not say which
control goes with which tumour. Getting that table would make the differential
analysis more powerful.

**A defect in how the model picks negative examples.** The code takes the first
50 kinases alphabetically rather than a random 50, so 87% of kinases are never
used as a negative example. We checked whether this biased the results and it
did not, so it is a flaw to document rather than a reason to redo anything. Any
fix would mean rerunning every experiment.

**Chapters I and VI of the thesis** are on hold until October 31, because the
framing depends on the first decision above.
