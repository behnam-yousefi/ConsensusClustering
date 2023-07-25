## Robust Stratification Through Consensus Clustering (CC)
## Type 3: Meta-data CC (MCC): Different views or obseravtions of the same phenomenon

rm(list = ls())

setwd("~/Desktop/R_Root/ConsensusCluster/")

## Definitions
library(cluster)
source("R/cluster_generation_functions.R")
source("R/CC_functions.R")

## Data generation

set.seed(1)
data = multiview_clusters (n = c(100,100), hidden.dim = 2, observed.dim = c(2,2,2), 
                           sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))

X_observation = data[["observation"]]
X_hidden = data[["hidden"]]
class = data[["class"]]

plot(X_hidden[,1:2], pch = 20, col = class, cex = .5)
plot(X_observation[[3]][,1:2], pch = 20, col = class, cex = .5)


## MCC
Clusters = multiview_kmeans(X_observation, rep = 100, range.k = c(5,20), method = "random")
Clusters = multiview_pam(X_observation, rep = 100, range.k = c(5,20), method = "random")

Adj = coCluster_matrix(Clusters)
pheatmap::pheatmap(Adj)


CM = consensus_matrix(Adj, max.cluster = 5, resample.ratio = 0.8, max.itter = 50, clustering.method = "hclust")

Scores = CC_cluster_count(CM)

# Method 2
Clusters = multiview_pam(X_observation, rep = 1, range.k = c(2,10), method = "silhouette")
Adj = coCluster_matrix(Clusters)
pheatmap::pheatmap(1-Adj)

Clusters = multi_pam (1-Adj, rep = 1, range.k = c(2,5), is.distance = TRUE, method = "silhouette")
  

# Plots
RobScore = Scores[["RobScore"]]
plot(RobScore, pch = 20, type = "b", ylab = "robustness score", xlab = "number of clusters")
RobScore = Scores[["PAC"]]
plot(RobScore, pch = 20, type = "b", ylab = "PAC", xlab = "number of clusters")

Kopt = Scores[["Kopt_RobScore"]]
clusters = HirClustFromAdjMat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE, method = "ward.D")
# clusters = SpectClustFromAdjMat(CM[[Kopt]], k = Kopt, max.eig = Kopt, alpha = 1, adj.conv = FALSE)
plot(X_hidden[,1:2], pch = 20, col = clusters, cex = .5)
plot(X_observation[[1]][,c(1,2)], pch = 20, col = clusters, cex = .5)

