## Robust Stratification Through Consensus Clustering (CC)
## Type 1: Permuted CC (PCC): Different permutations of samples
## Implementations of the PCC and compare to baselines

rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusCluster/")

## Definitions
library(cluster)
source("R/cluster_generation_functions.R")
source("R/CC_functions.R")

## Data generation
dim  = 2
data = gaussian_clusters(n = c(40,40,40), dim = dim, sd.max = .5, sd.noise = 0, r.range = c(1,2))
X = data$X
class = data$class
plot(X[,1:2], pch = 20, col = class, cex = .5)


## PCC
dists = as.matrix(dist(X))
Adj = (max(dists) - dists)/max(dists)

X = X / apply(X, 1, function(x){sqrt(sum(x^2))})
Cos = X %*% t(X)
Adj = Cos

CM = consensus_matrix(Adj, max.cluster = 5, resample.ratio = 0.8, max.itter = 50, clustering.method = "hclust")

Scores = CC_cluster_count(CM)
RobScore = Scores[["PAC"]]
plot(RobScore, pch = 20, type = "b", ylab = "robustness score", xlab = "number of clusters")

Kopt = Scores[["Kopt_LogitScore"]]
pheatmap::pheatmap(CM[[2]])

clusters = HirClustFromAdjMat(CM[[Kopt]], k = Kopt, alpha = 1, adj.conv = FALSE)
clusters = SpectClustFromAdjMat(Adj, k = Kopt, max.eig = Kopt, alpha = 1, adj.conv = TRUE)
clusters = PamClustFromAdjMat (CM[[Kopt]], k = 2, alpha = 1, adj.conv = TRUE)

plot(X[,1:2], pch = 20, col = clusters, cex = .5)

coClMat = coCluster_matrix(cbind(clusters, data[,3]))
ACC = sum(coClMat==0) / (2*N)
print(ACC)


rownames(X) = paste0("Y",1:nrow(X))
res = M3C(t(X), maxK=5, clusteralg = "hc", objective = "PAC")
plot(X[,1:2], col=res$assignments)


## Other methods
Tree = hclust(as.dist(dists), method="complete")
Sil = rep(0,5)
for (k in 2:5){
  clusters = cutree(Tree, k=k)
  Sil[k] = mean(silhouette(clusters, as.dist(dists))[,3])
}
plot(Sil)
clusters = cutree(Tree, k=3)
plot(data[,1:2], pch = 20, col = clusters, cex = .5)
coClMat = coCluster_matrix(cbind(clusters, data[,3]))
ACC = sum(coClMat==0) / (N^2)
print(ACC)

Sil = rep(0,5)
for (k in 2:5){
  clusters = kmeans(X[,1:2], k)[["cluster"]]
  Sil[k] = mean(silhouette(clusters, as.dist(dists))[,3])
}
plot(Sil)
clusters = kmeans(X[,1:2], 3)[["cluster"]]
plot(data[,1:2], pch = 20, col = clusters, cex = .5)
coClMat = coCluster_matrix(cbind(clusters, data[,3]))
ACC = sum(coClMat==0) / (N*2)
print(ACC)

# A = Adj
# Th = .95
# A[A<Th] = 0
# A[A>=Th] = 1
clusters = SpectClustFromAdjMat(Adj, k = Kopt, max.eig = Kopt, alpha = 1, adj.conv = TRUE)
plot(data[,1:2], pch = 20, col = clusters, cex = .5)
coClMat = coCluster_matrix(cbind(clusters, data[,3]))
ACC = sum(coClMat==0) / (N*2)
print(ACC)

## M3C
library(corpcor)
library(ggplot2)
X = clustersim(n=100, n2 = 2, r = 5, K = 3, alpha = .6, wobble = .5)[["simulated_data"]]
X = t(X)
plot(X[,1:2])

