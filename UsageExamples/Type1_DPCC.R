## Robust Stratification Through Consensus Clustering (CC)
## Type 1: Permuted CC (PCC): Different permutations of samples
## Implementations of the PCC and compare to baselines
rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusCluster/")
library(ConsensusCluster)

## Data generation
dim  = 2
data = gaussian_clusters(n = c(40,40,40), dim = dim, sd.max = .5, sd.noise = 0, r.range = c(1,2))
X = data$X
class = data$class
plot(X[,1:2], pch = 20, col = class, cex = .5)

## Calculate the adjacency matrix
Adj = adj_mat(X, method = "euclidian")

## Apply DPCC

# 1) Generation mechanism and calculate consensus matrix
CM = consensus_matrix(Adj, max.cluster = 3, resample.ratio = 0.8, max.itter = 10, clustering.method = "hclust")

# 2) Calculate stability score and count the optimum number of clusters
Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, pch = 20, type = "b", ylab = "robustness score", xlab = "number of clusters")

# 3) Perform the final clustering
Kopt = Scores[["Kopt_LogitScore"]]
pheatmap::pheatmap(CM[[3]])

clusters = hir_clust_from_adj_mat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
clusters = spect_clust_from_adj_mat(Adj, k = Kopt, max.eig = Kopt, alpha = 1, adj.conv = TRUE)
clusters = pam_clust_from_adj_mat(CM[[Kopt]], k = 2, alpha = 1, adj.conv = TRUE)

plot(X[,1:2], pch = 20, col = clusters, cex = .5)

