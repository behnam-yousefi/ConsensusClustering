## Consensus clustering functions
## version that works for multiview data in cases where not all the data is available for all the datasets.
## Author: Behnam Yousefi

## Basic Functions -------------------------------------------------------------
#' Logit function
#'
#' @param x numerical scaler input
#' @value Logit(x) = log(1*x/(1-x))
#' @examples
#' y = Logit(0.5)
#'
Logit = function(x)
  return(log(1*x/(1-x)))


#' Build connectivity matrix
#'
#' @param clusters a vector of clusterings. Zero elements mean that the sample was absent
#' during clustering
#'
#' @value Connectivity matrix
#'
#' @details
#'  Connectivity matrix (M) is a binary matrix N-by-N
#'  M[i,j] = 1 if sample i and j are in the same cluster
#'  ref: Monti et al. (2003) "Consensus Clustering: A Resampling-Based Method for
#'  Class Discovery and Visualization of Gene Expression Microarray Data", Machine Learning
#'
#' @examples
#' connectivity_mat([1,1,1,2,2,2])
#'
connectivity_matrix = function (clusters){
  Nsample = length(clusters)
  M = matrix(0,Nsample,Nsample)
  for (i in 1:Nsample)
    for (j in 1:Nsample)
      if (clusters[i]*clusters[j]>0)
        M[i,j] = ifelse(clusters[i]==clusters[j],1,0)
  return(M)
}

#' Build indicator matrix
#'
#' @param clusters a vector of clusterings. Zero elements mean that the sample was absent
#' during clustering
#'
#' @value Indicator matrix
#'
#' @details
#'  Indicator matrix (I) is a binary matrix N-by-N
#'  I[i,j] = 1 if sample i and j co-exist for clustering
#'  ref: Monti et al. (2003) "Consensus Clustering: A Resampling-Based Method for
#'  Class Discovery and Visualization of Gene Expression Microarray Data", Machine Learning
#'
indicator_matrix = function (clusters){
  Nsample = length(clusters)
  I = matrix(0,Nsample,Nsample)
  for (i in 1:Nsample)
    for (j in 1:Nsample)
      I[i,j] = ifelse(clusters[i]*clusters[j]>0,1,0)
  return(I)
}

#' Convert adjacency function to the affinity matrix
#'
#' @param adj.mat Adjacency matrix. The elements must be within [-1, 1].
#' @param alpha soft threshold value (see details).
#'
#' @details
#' adj = exp(-(1-adj)^2/(2*alpha^2))
#' ref: Luxburg (2007), "A tutorial on spectral clustering", Stat Comput
#' @examples
#' Adj_mat = rbind(c(0.0,0.9,0.0),
#'                 c(0.9,0.0,0.2),
#'                 c(0.0,0.2,0.0))
#' adj_conv(Adj_mat)
#'
#'
adj_conv = function(adj.mat, alpha=1){
  aff_mat = exp(-(1-adj.mat)^2/(2*alpha^2))
  return(aff_mat)
}

#' Hierarchical clustering from adjacency matrix
#'
#' @param adj.mat adjacency matrix
#' @param k number of clusters (default=2)
#' @param alpha soft threshold (considered if \code{adj.conv = TRUE}) (default=1)
#' @param adj.conv binary value to apply soft thresholding (default=TRUE)
#'
#' @value vector of clusters
#'
#' @details
#' apply PAM (k-medoids) clustering on the adjacency matrix
#'
#' @examples
#' Adj_mat = rbind(c(0.0,0.9,0.0),
#'                 c(0.9,0.0,0.2),
#'                 c(0.0,0.2,0.0))
#' hir_clust_from_adj_mat(Adj_mat)
#'
hir_clust_from_adj_mat = function (adj.mat, k = 2, alpha = 1, adj.conv = TRUE, method = "ward.D"){
  adj.mat[is.na(adj.mat)] = 0

  Aff = if(adj.conv)
    Aff = adj_conv(adj.mat, alpha)
  else
    Aff = adj.mat

  dists = as.dist(1-Aff)
  Tree = hclust(dists, method=method)
  clusters = cutree(Tree, k=k)

  return(clusters)
}

#' Spectral clustering from adjacency matrix
#'
#' @param adj.mat adjacency matrix
#' @param k number of clusters (default=2)
#' @param alpha soft threshold (considered if \code{adj.conv = TRUE}) (default=1)
#' @param adj.conv binary value to apply soft thresholding (default=TRUE)
#'
#' @value vector of clusters
#'
#' @details
#' apply PAM (k-medoids) clustering on the adjacency matrix
#'
#' @examples
#' Adj_mat = rbind(c(0.0,0.9,0.0),
#'                 c(0.9,0.0,0.2),
#'                 c(0.0,0.2,0.0))
#' hir_clust_from_adj_mat(Adj_mat)
#'
spect_clust_from_adj_mat = function (adj.mat, k = 2, max.eig = 10, alpha = 1, adj.conv = TRUE, do.plot = FALSE){
  adj.mat[is.na(adj.mat)] = 0

  Aff = if(adj.conv)
    Aff = adj_conv(adj.mat, alpha)
  else
    Aff = adj.mat

  graph = graph_from_adjacency_matrix(Aff, mode="undirected", weighted=TRUE)
  Lsym = as.matrix(laplacian_matrix(graph, normalized = TRUE))
  Lsym[is.na(Lsym)] = 0

  PCA = prcomp(Lsym)
  Lambda = PCA$sdev[ncol(Lsym):1]
  Eigenvectors = PCA$rotation[,ncol(Lsym):1]

  Eigenvectors = Eigenvectors[,1:max.eig]
  clusters = kmeans(Eigenvectors, k)[["cluster"]]

  if (do.plot)
    plot(Lambda[1:10], pch = 20, type = "b")
  return(clusters)
}

#' PAM (k-medoids) clustering from adjacency matrix
#'
#' @param adj.mat adjacency matrix
#' @param k number of clusters (default=2)
#' @param alpha soft threshold (considered if \code{adj.conv = TRUE}) (default=1)
#' @param adj.conv binary value to apply soft thresholding (default=TRUE)
#'
#' @value vector of clusters
#'
#' @details
#' apply PAM (k-medoids) clustering on the adjacency matrix
#'
#' @examples
#' Adj_mat = rbind(c(0.0,0.9,0.0),
#'                 c(0.9,0.0,0.2),
#'                 c(0.0,0.2,0.0))
#' pam_clust_from_adj_mat(Adj_mat)
#'
pam_clust_from_adj_mat = function (adj.mat, k = 2, alpha = 1, adj.conv = TRUE){
  adj.mat[is.na(adj.mat)] = 0

  Aff = if(adj.conv)
    Aff = adj_conv(adj.mat, alpha)
  else
    Aff = adj.mat

  dists = as.dist(1-Aff)
  clusters = pam(as.matrix(dists), k = k, diss = TRUE, cluster.only=TRUE)

  return(clusters)
}


#' Calculate the  Co-cluster matrix for a set of given set of clustering results.
#'
#' @param X clustering matrix of Nsamples x Nclusterings.
#' Zero elements are are considered as unclustered samples
#'
#' @value The normalized matrix of Co-cluster frequency of any pairs of samples (Nsamples x Nsamples)
#'
#' @examples
#' Clustering = cbind(c(1,1,1,2,2,2),
#'                    c(1,1,2,1,2,2))
#' coCluster_matrix(Clustering)
#'
coCluster_matrix = function(X){
  Nsample = nrow(X)
  Nmethod = ncol(X)

  M = matrix(0,Nsample,Nsample)
  I = M
  pb = txtProgressBar(min = 0, max = Nmethod, style = 3)

  for (cl in 1:Nmethod){
    M = M + connectivity_matrix(X[,cl])
    I = I + indicator_matrix(X[,cl])
    setTxtProgressBar(pb, cl)
  }
  coClusterMatrix = M/I
  close(pb)

  rownames(coClusterMatrix) = rownames(X)
  colnames(coClusterMatrix) = rownames(X)
  return(coClusterMatrix)
}
