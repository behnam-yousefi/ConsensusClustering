## Robust Stratification Through Consensus Clustering (CC)
## Comparison framework

### read SLURM_ARRAY_TASK_ID
task_id = Sys.getenv("SLURM_ARRAY_TASK_ID")
task_id = as.numeric(task_id)
# task_id: XY, X = 10*sd_max, Y = rep
sd_max = floor(task_id/10)/10        #index for outer loop
rep = task_id%%10                    #index for iner loop
k = 3

print("TaskID is")
print(task_id)

setwd("~/ConsensusClustering/")

## Definitions
library(cluster)
source("cluster_generation_functions.R")
source("CC_functions.R")
source("M3C_functions.R")

## Parameters
Kopt = rep(0, 10)
names(Kopt) = c("PCC_PAM_RS", "PCC_PAM_PAC", "PCC_PAM_DA", 
                   "ECC_PAM_RS", "ECC_PAM_PAC", "ECC_PAM_DA",
                   "M3C_PAM_Ent_RCSI", "M3C_PAM_Ent_Pval", "Sil", "K_real")

print(dim)
print(rep)
print(k)
  
set.seed(task_id)
## Data generation
# data = gaussian_mixture_clusters(n = c(300,300,200), dim = dim, sd.max = .5, sd.noise = .1, r.range = c(1,2), mixture.range = c(3,4))
data = gaussian_clusters(n = c(20,20,20), dim = 2, sd.max = sd_max, sd.noise = 0, r.range = c(1,2))
X = data$X
class = data$class
plot(X[,c(1,2)], pch = 20, col = class, cex = .5)
  
## PCC
dists = as.matrix(dist(X))
Adj = (max(dists) - dists)/max(dists)
CM = consensus_matrix(Adj, max.cluster = 5, resample.ratio = 0.7, max.itter = 50, clustering.method = "pam", adj.conv = TRUE)
Scores = CC_cluter_count(CM)
Kopt[1] = Scores[["Kopt_RobScore"]]
Kopt[2] = Scores[["Kopt_PAC"]]
Kopt[3] = Scores[["Kopt_deltaA"]]
  
## ECC
Clusters = multi_kmeans(X, rep = 200, range.k = c(5,15), method = "random")
Adj = coCluster_matrix(Clusters)
CM = consensus_matrix(Adj, max.cluster = 10, resample.ratio = 0.7, max.itter = 50, clustering.method = "pam", adj.conv = FALSE)
  
Scores = CC_cluter_count(CM)
Kopt[4] = Scores[["Kopt_RobScore"]]
Kopt[5] = Scores[["Kopt_PAC"]]
Kopt[6] = Scores[["Kopt_deltaA"]]
  
## M3C
rownames(X) = paste0("Y",1:nrow(X))
res = M3C(t(X), maxK=5, clusteralg = "pam", objective = "entropy")
Kopt[7] = which.max(res[["scores"]]$RCSI) + 1
Kopt[8] = which.max(res[["scores"]]$P_SCORE) + 1
if (max(res[["scores"]]$P_SCORE)<1)
  Kopt[8] = 1
  
## Average Silhouette width
Sil = rep(0,10)
for (i in 2:10){
  clusters = kmeans(X, i)[["cluster"]]
  Sil[i] = mean(silhouette(clusters, as.dist(dists))[,3])
}
Kopt[9] = which.max(Sil)
  
Kopt[10] = k
print(Kopt)

# g: gaussian, p: pam, XYZ: K sd_max rep
saveRDS(Kopt, file = paste0("Data/temp/g_pam",k,task_id, ".rds"))