## Robust Stratification Through Consensus Clustering (CC)
## Type 3: Meta-data CC (MDCC): Different views or observations of the same phenomenon

rm(list = ls())
library(ConsensusClustering)

## Data generation
data = multiview_clusters (n = c(40,40,40), hidden.dim = 2, observed.dim = c(2,2,2),
                           sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))

X_observation = data[["observation"]]
X_hidden = data[["hidden"]]
class = data[["class"]]

plot(X_hidden, pch = 20, col = class, cex = .5)
plot(X_observation[[1]], pch = 20, col = class, cex = .5)
plot(X_observation[[2]], pch = 20, col = class, cex = .5)
plot(X_observation[[3]], pch = 20, col = class, cex = .5)


## Apply MCC

# Method 1) MDCC based on DPCC: each data view is considered as a perturbed version of the hidden state.

Adj = list()
for (i in 1:length(X_observation)){
  Adj[[i]] = adj_mat(X_observation[[i]], method = "euclidian")
}

CM = multiview_consensus_matrix(Adj, max.cluster = 5, clustering.method = "hclust")

Scores = CC_cluster_count(CM)
K_opt = Scores[["Kopt_LogitScore"]]
print(K_opt)


# Method 2) MDCC based on MPCC: perform MPCC on each view and then aggregate.

Clusters = multiview_kmeans_gen(X_observation, rep = 100, range.k = c(5,20), method = "random")
Clusters = multiview_pam_gen(X_observation, rep = 100, range.k = c(5,20), method = "random")

Adj = coCluster_matrix(Clusters)
pheatmap::pheatmap(Adj)

CM = consensus_matrix(Adj, max.cluster = 5, resample.ratio = 0.8, max.itter = 50, clustering.method = "hclust")
Scores = CC_cluster_count(CM)
K_opt = Scores[["Kopt_LogitScore"]]
print(K_opt)

