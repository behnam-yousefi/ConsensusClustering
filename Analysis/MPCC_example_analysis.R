## Robust Stratification Through Consensus Clustering (CC)
## Figure 3

rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusClustering/")

## Definitions
library(cluster)
source("R/random_data_generation.R")
source("R/CC_functions.R")

## Data generation
data = generateSpiralData(500)
plot(data[,1:2], pch = 20, col = data[,3], cex = .5)

# add point (5,0) to the data set
X = data[,1:2]
o = c(5,-.5)
X = rbind(o, X)
plot(X, pch = 20, col = ifelse(c(1,data[,3])==1, "red", "blue"), cex = .5)

# Figure A
clusters = kmeans(X, 2)[["cluster"]]
pdf(file = "Figure3A.pdf", 7,7)
plot(X, pch = 20, cex = .7, xlim = c(-5,10), ylim = c(-5,10),
     col = ifelse(clusters==1, "red", "blue"))
dev.off()

# Plot Euclidian distance from point O
O = cbind(rep(o[1], nrow(X)), rep(o[2], nrow(X)))
EucDist_from_O = apply(X-O , 1, function(x){return(sum(x^2))})
Col = colorRamps::matlab.like(10)
plot(X, pch = 20, cex = .6, xlim = c(-5,10), ylim = c(-5,10),
     col = Col[findInterval(1/EucDist_from_O, vec = seq(0,1,length.out=10), all.inside = TRUE)])

# Figure B
## ECC
Clusters = multi_kmeans(X, rep = 500, range.k = c(10,100), method = "random")
sim = coCluster_matrix(Clusters)
# pheatmap::pheatmap(sim)

d = as.dist(1-sim)
Tree = hclust(d, method = "single")
clusters = cutree(Tree, k=2)

pdf(file = "Figure3B.pdf", 7,7)
plot(X, pch = 20, cex = .7, xlim = c(-5,10), ylim = c(-5,10),
     col = ifelse(clusters==1, "red", "blue"))
dev.off()

# Figure C
CMDist_from_O = sim[1758,]
CMDist_from_O = sim[1763,]
CMDist_from_O = sim[1758,]

pdf(file = "Figure3C.pdf", 7,7)

Col = colorRamps::matlab.like(100)
Col = viridisLite::rocket(100)[100:1]
Col = grDevices::heat.colors(100)[100:1]
Col = rgb(sequence(100)/100,sequence(100)/100,sequence(100)/100)[100:1]

plot(X, pch = 21, cex = .6, col = "black" , lwd=.2, xlim=c(-5,10), ylim=c(-5,10),
     bg = Col[findInterval((CMDist_from_O)^.5, vec = seq(0,1,length.out=100), all.inside = TRUE)])
dev.off()