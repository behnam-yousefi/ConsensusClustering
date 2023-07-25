## Robust Stratification Through Consensus Clustering (CC)
## Type 1: Permutated CC (PCC): Different permutations of samples
## Implementations of the PCC and compare to baselines

rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusCluster/")

## Definitions
library(cluster)
source("R/cluster_gen_functions.R")
source("R/CC_functions.R")

## Data generation
SDnoise = 1
n = c(50, 50, 30)
center = rbind(c(5, 5), c(1,1), c(5,-1))
sigma = list(matrix(c(3, 1, 1, 2), nrow = 2),
             matrix(c(1, 0, 0, 1), nrow = 2),
             matrix(c(1, 0, 0, 1), nrow = 2))

data = generateGaussianClusters(n = n, center = center, sigma =  sigma)
plot(data[,1:2], pch = 20, col = data[,3], cex = .5)

dim = ncol(center)
N = sum(n)
data[,1:dim] = data[,1:dim] + matrix(rnorm(N*dim, 0, SDnoise), N)
plot(data[,1:2], pch = 20, col = data[,3], cex = .5)

## PCC
X = data[,1:dim]

Clusters = multi_kmeans(X, rep = 200, range.k = c(40,50), method = "random")
Adj = coCluster_matrix(Clusters)
pheatmap::pheatmap(Adj)

CM = consensus_matrix(Adj, max.cluster = 4, resample.ratio = 0.7, max.itter = 50, clustering.method = "hclust",
                      is.similarity = TRUE, no.cores = 1, adj.conv = FALSE)

Scores = CC_cluster_count(CM)
RobScore = Scores[["RobScore"]]
plot(RobScore, pch = 20, type = "b", ylab = "robustness score", xlab = "number of clusters")
pheatmap::pheatmap(CM[[3]])

Kopt = Scores[["Kopt_RobScore"]]

clusters = HirClustFromAdjMat(CM[[2]], k = Kopt, alpha = 1, adj.conv = FALSE, method = "single")
clusters = SpectClustFromAdjMat(Adj, k = Kopt, max.eig = Kopt, alpha = 1, adj.conv = FALSE)
plot(X[,1:2], pch = 20, col = clusters, cex = .5)

coClMat = coCluster_matrix(cbind(clusters, data[,3]))
ACC = sum(coClMat==0) / (2*N)
print(ACC)
