## Robust Stratification Through Consensus Clustering (CC)
## Type 1: Data perturbation CC (DPCC)

rm(list = ls())
library(ConsensusClustering)

## Data simulation
dim  = 2
data = gaussian_clusters(n = c(40,40,40), dim = dim, sd.max = .5, sd.noise = 0, r.range = c(1,2))
X = data$X
class = data$class
plot(X, pch = 20, col = class, cex = .5)

## Calculate the adjacency matrix (as the algorithm requires a similarity matrix)
Adj = adj_mat(X, method = "euclidian")

## Apply DPCC

# 1) Generation mechanism and calculate consensus matrix
CM = consensus_matrix(Adj, max.cluster = 5, resample.ratio = 0.7, max.itter = 50, clustering.method = "hclust")

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

