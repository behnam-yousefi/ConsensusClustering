## Consensus clustering with network aggregation
## By: Behnam Yousefi

rm(list = ls())
setwd("~/Desktop/R_Root/ConsensusClustering/")
source("R/CC_functions.R")
source("R/ActivityPlot.R")
library(tmod)
library(ggplot2)

### Find highly variable genes across all the train datasets
## 1-Find data set specific highly variable gebes
data_root = "~/Desktop/R_Root/SLEmap/Data/Curated_Data/"
gseID_Train = c("GSE65391", "GSE49454", "GSE45291", "GSE22098", "GSE39088")
GeneNamesToKeep = list()
for (gseID in gseID_Train) {
  pheno = readRDS(paste0(data_root,gseID,".RDS"))[["phenotype"]][["curated"]]
  expr = readRDS(paste0(data_root,gseID,".RDS"))[["expr"]]
  print(nrow(expr))
  expr = expr[!duplicated(rownames(expr)),]
  print(nrow(expr))
  # tokeep = order(apply(expr, 1, sd), decreasing=TRUE)[1:5000]
  tokeep = order(apply(expr[,pheno$disease %in% c("SLE","ASLE","PSLE")], 1, sd), decreasing=TRUE)[1:5000]
  # tokeep = order(apply(expr, 1, function(x){
  #   p_val = wilcox.test(x[pheno$disease %in% c("SLE","ASLE","PSLE")], x[!pheno$disease %in% c("SLE","ASLE","PSLE")])$p.value
  #   return(-log(p_val))}), decreasing=TRUE)[1:5000]
  tokeepNames = rownames(expr[tokeep,])

  tokeepNames = tokeepNames[!is.na(tokeepNames)]
  tokeepNames = tokeepNames[! tokeepNames %in% ""]
  tokeepNames = tokeepNames[! tokeepNames %in% " "]

  GeneNamesToKeep[[gseID]] = tokeepNames
}

## 2-Find intersected highly variable genes across the train data set
UnionGenes = GeneNamesToKeep[[gseID_Train[1]]]
for (i in 2:length(gseID_Train))
  UnionGenes = union(UnionGenes, GeneNamesToKeep[[gseID_Train[i]]])

length(UnionGenes)

UnionGenes = data.frame(genes = UnionGenes, rep = rep(0,length(length(UnionGenes))))
for (i in 1:length(gseID_Train)){
  genes_i = GeneNamesToKeep[[gseID_Train[i]]]
  UnionGenes$rep[UnionGenes$genes %in% genes_i] = UnionGenes$rep[UnionGenes$genes %in% genes_i] + 1
}

SelectedGenes = UnionGenes$genes[UnionGenes$rep>=3]
length(SelectedGenes)

## 3-Build individual gene co-expression networks with the intersected genes
AdjMat = list()
Nsamples = list()
for (gseID in gseID_Train) {
  pheno = readRDS(paste0(data_root,gseID,".RDS"))[["phenotype"]][["curated"]]
  expr = readRDS(paste0(data_root,gseID,".RDS"))[["expr"]]
  expr = expr[!duplicated(rownames(expr)),]
  expr = data.frame(expr)
  #expr = expr[GeneNamesToKeep[[gseID]],]
  expr = expr[GeneNamesToKeep[[gseID]],pheno$disease %in% c("SLE","ASLE","PSLE")]
  expr = t(scale(t(expr)))
  AdjMat[[gseID]] = cor(t(expr))
  AdjMat[[gseID]][is.na(AdjMat[[gseID]])] = 0
  Nsamples[[gseID]] = ncol(expr)
}

## 4-Consensus Clustering

# CM = multiview_consensus_matrix(AdjMat, max.cluster = 4, sample.set = SelectedGenes, clustering.method = "pam")
# saveRDS(CM, file = "Data/CM_MCC_SLE_final2.rds")
CM = readRDS("Data/CM_MCC_SLE_final.rds")
Scores = CC_cluster_count(CM)

K_opt = Scores[["Kopt_LogitScore"]]
print(K_opt)

# pheatmap::pheatmap(CM[[K_opt]])


clusters_list = list()
clusters = PamClustFromAdjMat (CM[[K_opt]], k = K_opt, alpha = 1, adj.conv = TRUE)
names(clusters) = SelectedGenes
clusters_list[[1]] = clusters
for (i in 1:5){
  clusters = PamClustFromAdjMat (AdjMat[[i]], k = K_opt, alpha = 1, adj.conv = TRUE)
  clusters = clusters[SelectedGenes]
  table(clusters)
  clusters_list[[i+1]] = clusters
}

### Figure 6
PW_all = list()
ID_all = c()
ID_all_fullname = c()
for (j in 1:6){
  clusters = clusters_list[[j]]
  PW_all[[j]] = list()
  for (i in 1:4){
    genesK = names(clusters[clusters==i])
    pw = tmodUtest(genesK,qval=0.05, mset="all", useR=FALSE)
    PW_all[[j]][[i]] = pw

    ID_all_fullname = c(ID_all_fullname, paste0(pw$Title, " (", pw$ID,")"))
    ID_all = c(ID_all, pw$ID)
  }
}
ID_all_fullname = unique(ID_all_fullname)[-61]
ID_all = unique(ID_all)

j = 5
PW_mat = data.frame(matrix(1, length(ID_all), K_opt))
rownames(PW_mat) = ID_all
for (i in 1:K_opt){
  pw_i = PW_all[[j]][[i]]
  PW_mat[rownames(pw_i), i] = pw_i$adj.P.Val
}
max(-log10(PW_mat))


# Order = c(1, 4:7, 3, 9, 8, 11,18,32,33,67,  12,14,19,20,23,28,31,43,47,53,  44,13,16,24,29,34,35,37,39,   21,36,46,48,50,49,51,52,59,66,74,   54,55,56,58,       2,10,26,27,30,38,15,17,40,45,22,25,41,42,57,60,61:65,68:73)
# Order = c(1, 4:7, 3, 9, 8, 11,18,32,33,67,  12,14,19,20,23,28,31,  44,13,16,24,29,34,35,37,39,   21,36,46,48,50,49,51,52,59,66,74,   54,55,56,57,  10,26,27,30,38,15,17,40,45,22,25,41,42,58,60,61:65,68:73)
Order = c(46,48,50,49,51,52,59,66,  1, 4:7, 3, 9, 8, 54,55,56, 11,18,32,33,67,  12,14,19,20,23,28,31,53,  44,13,16,24,29,34,35,37,39,74, 21,36,    10,26,27,30,38,15,17,40,45,22,25,41,42,58,57,60,61:65,68:73)

ID_all_fullname_sorted = ID_all_fullname[Order]
PW_mat_sorted = PW_mat[Order, c(3,1,4,2)]

# pdf("Figure6f_new.pdf")
activity_plot(-log10(PW_mat_sorted), PW_mat_sorted, alpha.min = 0, shape="c", size=c(1,3), p.val.th = 0.05,
              title = gseID_Train[j-1], pw.names=ID_all_fullname_sorted, h.line=TRUE, v.line=FALSE)
# dev.off()

# Figure 6.5
# Do this for each dataset and save "N_uniqe_pw" and "N_dup_pw"
i = 6
{
clusters = clusters_list[[i]]
k = 4
PW = list()
for (i in 1:k){
  genesK = names(clusters[clusters==i])
  pw = tmodUtest(genesK,qval=0.05, mset="all", useR=FALSE)
  PW[[i]] = pw
}
names(PW) = c("1","2","3","4")
#tmodPanelPlot(PW, text.cex = .5, min.e = 0.1, max.e = 1, plot.cex = .5, grid = "at", pval.thr=.05)

all_pw = sapply(PW , function(x)(return(x[["ID"]][x[["adj.P.Val"]]<0.05]) ))
pw_set = c(all_pw[[1]], all_pw[[2]], all_pw[[3]], all_pw[[4]])

N_uniqe_pw = length(unique(pw_set))
N_dup_pw = sum(!duplicated(pw_set[duplicated(pw_set)]))
print(N_uniqe_pw)
print(N_dup_pw)
}
# next data

N1 = c(32, 44, 42, 36, 33 ,57)
N2 = c( 2,  1,  4,  5,  3 , 3)
N2*100/N1

data_bar = data.frame(N = c(N1,N2),
                      dataset = as.factor(rep(c(gseID_Train, "Consensus\nClsutering"),2)),
                      cls = as.factor(c(rep("No. biological processes",6), rep("No. duplicated biological processes",6))))

pdf("Figure6g.pdf")
ggplot(data_bar, aes(x = dataset, y = N, fill = cls)) +
  geom_bar(stat = "identity", position = "dodge", width = .5)+
  scale_x_discrete(guide = guide_axis(angle = 90)) +
  scale_fill_manual(values=c("plum", "sienna2")) + theme_classic() + ylim(0,60)
dev.off()
