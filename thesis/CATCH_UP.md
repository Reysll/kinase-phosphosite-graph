# Catch-up: where this project stands

Written to be explained out loud. Last updated 2026-10-27.

## (a) What we set out to do

Kinases are enzymes that attach phosphate groups to specific spots on other
proteins. Those spots are phosphosites, and which kinase hits which site is what
we want to predict. Most of what is known comes from a general-purpose database
that is not specific to any tissue or disease, and more than 95% of known
phosphosites have no kinase assigned at all.

Our idea was to make the prediction context-aware. We built a network of
proteins and phosphosites, with edges for pathway membership, protein
interaction and co-evolution, then added liver cancer data by connecting
phosphosites whose levels move together across patients. The hypothesis was that
liver-specific information would improve predictions for liver cancer, and that
we could hand a short list of kinase and substrate pairs to a wet lab.

## (b) The six things we found

**1. The liver data does not actually add information.** We rebuilt the liver
network with random edges that keep each node's number of connections but
destroy the biology, and the fake network gives the same answer: 96% of
predictions come out with an identical rank. It also reproduced, rank for rank,
the one biologically exciting result we had, a cluster of necroptosis kinases
that looked more relevant in the cancer networks.

**2. Changing a random seed moves the result more than changing the biology
does.** Rerunning the same configuration under different seeds gave a
seed-to-seed difference of 0.373, on a scale where generic against liver is
0.249. Most of what we reported as a biological effect sits inside the noise of
rerunning the same experiment.

**3. The model barely looks at the phosphosite.** Of all the variation in a
kinase's predicted rank, 99.97% is explained by which kinase it is rather than
which site it is being scored against. The model produces one fixed ranking of
kinases and repeats it at every site, which is not site-specific prediction.

**4. A baseline with no model in it beats everything we built.** Ranking kinases
by how many substrates they already have on record, ignoring the network
entirely, gives a median rank of 25 against 74 for our best model. The one real
exception matters: for poorly studied kinases, the ones we most want to say
something about, our models win by a wide margin.

**5. We found a real defect, fixed it, and the fix made the core problem
worse.** The code picked its negative training examples alphabetically, so 87% of
kinases were never used as a negative, and a proper random draw did shift the
global ranking, by 5 median positions. But site dependence went the wrong way:
the 99.97% above rose to 99.998%, and the typical kinase now moves only 0.16
rank positions across all 3,228 sites.

**6. The differential analysis works, and it is the part we can hand over.**
Comparing phosphosite levels directly between tumour and adjacent healthy tissue
gives 659 sites that change substantially and reliably. It uses only the
measurements, not the network or the model, so none of the problems above apply
to it.

## (c) Why that happens, the leading explanation

We now have a measured answer, and it is not the one we expected.

**The model is sorting kinases by how many connections they have in the graph,
and handing back that same order at every site.** The correlation between a
kinase's connection count and its average predicted rank is **-0.9995**. The
substrate count from finding 4 manages only -0.44, so this is a different and
far stronger effect. The clearest case: PIK3R1 is ranked first at essentially
all 3,228 sites while having **two** known substrates, because it is the most
connected node in the graph. CDK1, with 693 substrates, comes third.

Connection count is in turn mostly pathway annotation: **95.5%** of the edges
touching a kinase are MSigDB pathway co-membership. So the model is really
ranking kinases by how many pathways they are annotated to.

We traced this to a bug, fixed it, reran everything, and the ranking did not
change. The embedding discarded one uninformative vector where it should have
discarded all of them, which left 16 of its 32 dimensions as the connection
count rescaled. Those are now gone, and the ranking is still connection count,
at -0.98 instead of -0.9995, with the median rank unmoved at 79 and 74. The
reason is measured: the count is not only in the vectors we dropped, it is in
the length of every node's vector, which is a property of this kind of embedding
rather than of our mistake. Making every vector the same length is the next
thing to try and is untested.

So the bug was real, the diagnosis of the ranking was wrong, and the honest
version is that connection count predicts our ranking at -0.98 and we do not yet
know how to stop it. Note the pattern with finding 5: two separate real defects,
both fixed properly, both fixes making the site-blindness worse. Detail in
`thesis/EMBEDDING_DEFECT.md`, section 9.

## (d) What the advisor asked for on October 15

| Request | Status |
|---|---|
| Journal paper, not a conference paper | Adopted as the direction |
| Significant liver-cancer substrates for the collaborator | Delivered, 659 sites |
| Comparative within-category ranking instead of accuracy | Done, result is negative |
| Re-do the preprocessing after the DEPOD count looked wrong | Audited, he was right |
| Explain the correlation sentence, meet step by step | Written up, meeting pending |

Three need nuance. The **category ranking** was already built back in June, so
we could test it the same day, and it does not separate the real network from
the random one: they agree at a correlation of 0.9944. On **preprocessing** he
was right, 115 phosphatases rather than 250, but that layer is 0.03% of the
network; the audit found something larger, that 90.6% of the network is pathway
co-membership, an annotation rather than a measurement, which is the likeliest
reason a random network does as well as the real one. On **correlation**, we do
not use all 18 samples: each pair uses only the samples where both sites were
detected, most often the minimum of 6, and only 8% of edges use all 18.

## (e) What is undecided

**What the main result of the paper is.** The ranking work cannot distinguish
itself from a random network, which is a legitimate negative finding but not the
paper we set out to write. Either the network and the features change, or the
differential phosphosite work becomes the main contribution and the ranking work
is reported as a cautionary result.

**What to try next, now that the obvious fix has failed.** Two candidates, both
touching project code, both waiting on you. One, make every node's vector the
same length, which removes the channel the first fix left behind. Two, tighten
the threshold on the correlation edges, which on measurement is what stops the
sites looking identical to the model, and which the failed fix does not touch. I
would do the second first, because it is the one whose target has not already
resisted a repair. Option D, the relational-only feature idea, is probably not
the lever either, since it cannot recover site information that is not in the
embedding to start with.

**Whether to change the network and the sampler.** The pathway edges at 90.6%
are now implicated directly. Making the random negative sampler the default is
still the better method, but it means rerunning the other two embedding methods
for consistency, roughly a day of compute.

**Where the liver data came from.** The file holds supplementary tables S4 and
S11 of an unidentified hepatocellular carcinoma paper. Until that is resolved we
cannot cite the dataset or get the patient-level information, including which
control belongs to which tumour, which would make the differential analysis more
powerful.

**Chapters I and VI** are on hold until October 31, because the framing depends
on the first decision above.
