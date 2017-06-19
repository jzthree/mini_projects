# scripts
This repository collects my mini-scripts which usually address some mini problems I don't know of an existing satisfying solution for. These are mostly experimenta. Any contribution / discussion is very welcomed (if anybody is reading this haha).


### order.tsp.R: 
Heatmap ordering is usually done with hierarchical clustering. Hierarchical clustering are not designed to give an optimal ordering so it is often tricky to get a good visualization. This script find an 'optimal' order the by minimizing the sum of distances between neighbors in the heatmap. This problem is equivalent to a Hamiltonian path problem which reduce to a special case of traveling salseman problem. So this script calls a binary of a TSP solver Concorde to find the ordering. It is usually very fast and give good vis for me. It is especially suitable for data with sequential structure. 


### poissonLikDist.py: 
many distance metric I previously knew are not particularly good for count data that are not normalized, especially when the count number is low. Specifically RNA-seq with low depth is the scenario that I use this for. Distance metic should indicate the uncertainty of the data so lower counts should generally result in higher distances and this is an attempt.
