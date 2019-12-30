# Mini-project scripts
This repository collects my scripts each addresses a mini problem for which I don't know any pre-existing satisfying solution. All of these are tested but some are very experimental so use at your own risk. Any contribution / discussion is very welcomed.

### ndea.R:
This implements a two sample weighted t-test where weighted averages rather than averages are compared. Weighted t-tests are widely considered before but it is worth noting that they generally consider the case where weights indicate equivalent sample sizes or consider weights as random variables. Here each weight is considered a fixed, *non-random*, property tied to each sample. When the weights are derived from a network (such as using network connectivity to a gene or gene set of interest) we call it NDEA (Network differential enrichment analysis).

This is published in Zhou J, Park CY, Theesfeld CL, Wong AK, Yuan Y, Scheckel C, Fak JJ, Funk J, Yao K, Tajima Y, Packer A, Darnell RB, Troyanskaya OG. (2019). Whole-genome deep-learning analysis identifies contribution of noncoding mutations to autism risk. Nature Genetics.


### order.tsp.R: 
Heatmap ordering is usually done with hierarchical clustering. Hierarchical clustering is not designed to give an optimal ordering so it is often tricky to get a good visualization. This script find an 'optimal' order the by minimizing the sum of distances between neighbors in the heatmap. This problem is equivalent to a Hamiltonian path problem which reduce to a special case of traveling salseman problem. So this script calls a binary of a TSP solver Concorde to find the ordering. It is usually very fast and give good vis for me. It is especially suitable for data with sequential structure. 


### poissonLikDist.py: 
Most distance functions are not particularly designed for count data that are not normalized (the sum of each sample varies but we only care about proportions). Distance metic should indicate the uncertainty of the data so lower counts should generally result in higher distances and this is an attempt. Specifically RNA-seq with low depth is the scenario that I use this for.


### weighted.cortest.R:

Sometime it is useful to weight samples differently when computing correlation, but surprisingly I have not found an implementation that calculate significance for the weighted correlation in the case that the weights are known and fixed (i.e. not drawn from a distribution) (note this is different from the weights R package). The trick is simply replacing sample size with the effective sample size derived from the weights.
