## Robust Stratification Through Consensus Clustering (CC)
## Type 2: Method perturbation CC (MPCC)

rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusClustering/")
library(ConsensusClustering)

## Data simulation
dim  = 2
data = gaussian_clusters(n = c(40,40,40), dim = dim, sd.max = .5, sd.noise = 0, r.range = c(1,2))
X = data$X
class = data$class
plot(X, pch = 20, col = class, cex = .5)

## Apply MPCC

# 1) Generation mechanism: calculate different clustering
Clusters = multi_kmeans_gen(X, rep = 200, range.k = c(2,10), method = "random")
# or
Clusters = multi_pam_gen(X, rep = 200, range.k = c(2,10), method = "random")

# 2) Consensus mechanism
# 2.1) Consensus mechanism using co-association matrix
Adj = coCluster_matrix(Clusters)
pheatmap::pheatmap(Adj)

# We can then use DPCC on the co-association matrix and count the number of clusters
CM = consensus_matrix(Adj, max.cluster = 5, resample.ratio = 0.7, max.itter = 50, clustering.method = "hclust")

# Calculate stability score and count the optimum number of clusters
Scores = CC_cluster_count(CM)
RobScore = Scores[["LogitScore"]]
plot(RobScore, pch = 20, type = "b", ylab = "robustness score", xlab = "number of clusters")

# 3) Perform the final clustering
Kopt = Scores[["Kopt_LogitScore"]]
pheatmap::pheatmap(CM[[Kopt]])

