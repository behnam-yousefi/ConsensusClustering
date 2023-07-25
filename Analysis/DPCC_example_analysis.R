## Figure 2
## Stability score explanation

rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusCluster/")

## Definitions
library(cluster)
source("Other/M3C_functions.R")
source("R/cluster_generation_functions.R")
source("R/CC_functions.R")

# Figure 2A: data generation
dim  = 2
set.seed(400)
data = gaussian_clusters(n = c(100,100,100,100), dim = dim, sd.max = .4, sd.noise = 0, r.range = c(1,2))
X = data$X
class = data$class

pdf("Figures/Figure_2A.pdf", 5,5)
plot(X[,c(1,2)], pch = 20, col = class, cex = 1)
plot(X[,c(1,2)], pch = 20, cex = 1)
dev.off()

## Figure 2B-F: CM
dists = as.matrix(dist(X))
Adj = (max(dists) - dists)/max(dists)
CM = consensus_matrix(Adj, max.cluster = 6, resample.ratio = 0.7, max.itter = 50, 
                      clustering.method = "pam", adj.conv = TRUE)

# saveRDS(CM, "Data/Fig_2_CM_seed400.rds")
CM = readRDS("Data/Fig_2_CM_seed400.rds")

jpeg("Figures/Figure_2B.jpeg")
pheatmap::pheatmap(CM[[2]], show_rownames=F, show_colnames=F, 
                   treeheight_row = 0, treeheight_col = 0,
                   clustering_method = "ward.D")
dev.off()
jpeg("Figures/Figure_2C.jpeg")
pheatmap::pheatmap(CM[[3]], show_rownames=F, show_colnames=F, 
                   treeheight_row = 0, treeheight_col = 0,
                   clustering_method = "ward.D")
dev.off()
jpeg("Figures/Figure_2D.jpeg")
pheatmap::pheatmap(CM[[4]], show_rownames=F, show_colnames=F,
                   treeheight_row = 0, treeheight_col = 0,
                   clustering_method = "ward.D")
dev.off()
jpeg("Figures/Figure_2E.jpeg")
pheatmap::pheatmap(CM[[5]], show_rownames=F, show_colnames=F, 
                   treeheight_row = 0, treeheight_col = 0,
                   clustering_method = "ward.D")
dev.off()
jpeg("Figures/Figure_2F.jpeg")
pheatmap::pheatmap(CM[[6]], show_rownames=F, show_colnames=F, 
                   treeheight_row = 0, treeheight_col = 0,
                   clustering_method = "ward.D")
dev.off()

## Figure 2G: cumulative distribution function
# modify "CC_cluster_count" for the Figure in paper
# K -> c(2,3,5,6,4)
# + ConsDistr = c(0, ConsDistr, 1)
# + x = c(0, x, 1)
# + y = c(0, y, 1)
# - text(k/10, plt[10*k], labels = k)
# - text(k/10, ConsDistr[10*k], labels = k)

pdf("Figures/Figure_2G_4.pdf", 5,5)
Scores = CC_cluster_count(CM, plot.cdf = TRUE, plot.logit = FALSE)
dev.off()

## Figure 2L: logit estimation of cumulative distribution function
pdf("Figures/Figure_2L_4.pdf", 5,5)
Scores = CC_cluster_count(CM, plot.cdf = FALSE, plot.logit = TRUE) 
dev.off()

## Figure 2H-K: robustness scores
pdf("Figures/Figure_2H.pdf",5,5)
RobScore = Scores[["deltaA"]]
RobScore[1] = NA
plot(RobScore, pch = 20, cex = 1.5, lwd = 1.5, type = "b", ylab = "deltaA score", xlab = "number of clusters")
dev.off()

pdf("Figures/Figure_2I.pdf",5,5)
RobScore = Scores[["PAC"]]
RobScore[1] = NA
plot(RobScore, pch = 20, cex = 1.5, lwd = 1.5, type = "b", ylab = "PAC score", xlab = "number of clusters")
dev.off()

pdf("Figures/Figure_2K.pdf",5,5)
RobScore = Scores[["LogitScore"]]
RobScore[1] = NA
plot(RobScore, pch = 20, cex = 1.5, lwd = 1.5, type = "b", ylab = "Logit score", xlab = "number of clusters")
dev.off()

# M3C
rownames(X) = paste0("Y",1:nrow(X))
res = M3C(t(X), maxK = 6, clusteralg = "pam", objective = "entropy")
saveRDS(res, file = "Data/Fig_2_M3C_seed400.rds")

pdf("Figures/Figure_2J.pdf",5,5)
RobScore = c(NA, res[["scores"]]$RCSI)
plot(RobScore, pch = 20, cex = 1.5, lwd = 1.5, type = "b", ylab = "RCSI score", xlab = "number of clusters")
dev.off()

