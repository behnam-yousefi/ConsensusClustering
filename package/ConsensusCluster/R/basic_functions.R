## Consensus clustering functions
## version that works for multiview data in cases where not all the data is available for all the datasets.
## Author: Behnam Yousefi

#' Covert data matrix to adjacency matrix
#'
#' @param X a matrix of samples by features.
#' @param method method for distance calculation:
#' \code{"euclidian"}, \code{"cosine"}, \code{"maximum"}, \code{"manhattan"},
#' \code{"canberra"}, \code{"binary"}, \code{"minkowski"},
#'
#' @return
#' calculated adjacency matrix from the data matrix using the specified methods
#'
#' @examples
#' X = gaussian_clusters()$X
#' Adj = adj_mat(X, method = "euclidian")
#'
adj_mat = function(X, method = "euclidian"){

  if (method == "cosine"){
    X = X / apply(X, 1, function(x){sqrt(sum(x^2))})
    AdjMat = X %*% t(X)
  }else{
    dists = as.matrix(stats::dist(X, method))
    AdjMat = (max(dists) - dists)/max(dists)
  }
  return(AdjMat)
}

#' Logit function
#'
#' @param x numerical scaler input
#' @return Logit(x) = log(1*x/(1-x))
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
#' @return Connectivity matrix
#'
#' @details
#'  Connectivity matrix (M) is a binary matrix N-by-N
#'  M[i,j] = 1 if sample i and j are in the same cluster
#'  ref: Monti et al. (2003) "Consensus Clustering: A Resampling-Based Method for
#'  Class Discovery and Visualization of Gene Expression Microarray Data", Machine Learning
#'
#' @examples
#' con_mat = connectivity_matrix(c(1,1,1,2,2,2))
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
#' @return Indicator matrix
#'
#' @details
#'  Indicator matrix (I) is a binary matrix N-by-N
#'  I[i,j] = 1 if sample i and j co-exist for clustering
#'  ref: Monti et al. (2003) "Consensus Clustering: A Resampling-Based Method for
#'  Class Discovery and Visualization of Gene Expression Microarray Data", Machine Learning
#'
#' @examples
#' ind_mat = indicator_matrix(c(1,1,1,0,0,1))
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
#' @param method distance meethod (default: /code{ward.D})
#'
#' @return vector of clusters
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

  dists = stats::as.dist(1-Aff)
  Tree = stats::hclust(dists, method=method)
  clusters = stats::cutree(Tree, k=k)

  return(clusters)
}

#' Spectral clustering from adjacency matrix
#'
#' @param adj.mat adjacency matrix
#' @param k number of clusters (default=2)
#' @param max.eig maximum number of eigenvectors in use (dafaut = 10).
#' @param alpha soft threshold (considered if \code{adj.conv = TRUE}) (default = 1)
#' @param adj.conv binary value to apply soft thresholding (default = \code{TRUE})
#' @param do.plot binary value to do plot (dafaut = \code{FALSE})
#'
#' @return vector of clusters
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

  graph = igraph::graph_from_adjacency_matrix(Aff, mode="undirected", weighted=TRUE)
  Lsym = as.matrix(igraph::laplacian_matrix(graph, normalized = TRUE))
  Lsym[is.na(Lsym)] = 0

  PCA = stats::prcomp(Lsym)
  Lambda = PCA$sdev[ncol(Lsym):1]
  Eigenvectors = PCA$rotation[,ncol(Lsym):1]

  Eigenvectors = Eigenvectors[,1:max.eig]
  clusters = stats::kmeans(Eigenvectors, k)[["cluster"]]

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
#' @return vector of clusters
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

  dists = stats::as.dist(1-Aff)
  clusters = cluster::pam(as.matrix(dists), k = k, diss = TRUE, cluster.only=TRUE)

  return(clusters)
}


#' Calculate the Co-cluster matrix for a given set of clustering results.
#'
#' @param X clustering matrix of Nsamples x Nclusterings.
#' Zero elements are are considered as unclustered samples
#' @param verbos binary value for verbosity (default = \code{TRUE})
#'
#' @return The normalized matrix of Co-cluster frequency of any pairs of samples (Nsamples x Nsamples)
#'
#' @details
#' Co-cluster matrix or consensus matrix (CM) is a method for consensus mechanism explaned in Monti et al. (2003).
#'
#' @examples
#' Clustering = cbind(c(1,1,1,2,2,2),
#'                    c(1,1,2,1,2,2))
#' coCluster_matrix(Clustering, verbos = FALSE)
#'
coCluster_matrix = function(X, verbos = TRUE){
  Nsample = nrow(X)
  Nmethod = ncol(X)

  M = matrix(0,Nsample,Nsample)
  I = M

  if (verbos)
    pb = utils::txtProgressBar(min = 0, max = Nmethod, style = 3)

  for (cl in 1:Nmethod){
    M = M + connectivity_matrix(X[,cl])
    I = I + indicator_matrix(X[,cl])

    if (verbos)
      utils::setTxtProgressBar(pb, cl)
  }
  coClusterMatrix = M/I

  if (verbos)
    close(pb)

  rownames(coClusterMatrix) = rownames(X)
  colnames(coClusterMatrix) = rownames(X)
  return(coClusterMatrix)
}

#' Similarity between different clusters
#'
#' @param x1 clustering vector 1
#' Zero elements are are considered as unclustered samples
#' @param x2 clustering vector 2
#' Zero elements are are considered as unclustered samples
#'
#' @return matrix of similarities between clustering labels
#'
#' @details
#' When performing performing several clustering, the cluster labels may no match with each other.
#' To find correspondences between clusters, the similarity between different labels need to be calculated.
#'
#' @examples
#' X = gaussian_clusters()$X
#' x1 = kmeans(X, 5)$cluster
#' x2 = kmeans(X, 5)$cluster
#' Sim = lebel_similarity(x1, x2)
#'
lebel_similarity = function(x1, x2){

  assertthat::assert_that(length(x1) == length(x2))

  clusters1 = sort(unique(x1))
  clusters1 = clusters1[clusters1>0]
  clusters2 = sort(unique(x2))
  clusters2 = clusters2[clusters2>0]

  assertthat::assert_that(length(clusters1) == length(clusters2))
  N_cluster = length(clusters1)

  Similarity = matrix(rep(0, N_cluster*N_cluster), N_cluster)
  rownames(Similarity) = clusters1
  colnames(Similarity) = clusters2

  for (i in clusters1){
    x1_bin = x1==i
    for (j in clusters2){
      x2_bin = x2==j
      sim = sum(x1_bin & x2_bin) / sqrt(sum(x1_bin) * sum(x2_bin))
      Similarity[i,j] = sim
    }
  }
  return(Similarity)
}

#' Relabeling clusters based on cluster similarities
#'
#' @param x1 clustering vector 1
#' Zero elements are are considered as unclustered samples
#' @param x2 clustering vector 2
#' Zero elements are are considered as unclustered samples
#'
#' @return dataframe of relabeled clusters
#'
#' @details
#' When performing performing several clustering, the cluster labels may no match with each other.
#' To perform maximum voting, the clustering need to be relabels based on label similarities.
#'
#' @examples
#' X = gaussian_clusters()$X
#' x1 = kmeans(X, 5)$cluster
#' x2 = kmeans(X, 5)$cluster
#' clusters = cluster_relabel(x1, x2)
#'
cluster_relabel = function(x1, x2){

  assertthat::assert_that(length(x1) == length(x2))

  ind2row_col = function(ind, nrow){
    col_index = ceiling(ind / nrow)
    row_index = ind %% nrow
    row_index = ifelse(row_index == 0, nrow, row_index)
    return(c(row_index, col_index))
  }

  Similarity = lebel_similarity(x1, x2)
  N_cluster = nrow(Similarity)

  index = order(Similarity, decreasing = TRUE)
  Map = c()
  for (i in index){
    map = ind2row_col(i, N_cluster)
    if (!(map[1] %in% Map[,1]) & !(map[2] %in% Map[,2]))
      Map = rbind(Map, c(map[1], map[2]))
  }

  Clusters = data.frame(x1 = x1, x2_old = x2, x2_new = rep(0, length(x2)))
  for (i in 1:N_cluster)
    Clusters$x2_new[Clusters$x2_old == Map[i,2]] = Map[i,1]

  return(Clusters)
}


#' Consensus mechanism based on majority voting
#'
#' @param X clustering matrix of Nsamples x Nclusterings.
#' Zero elements are are considered as unclustered samples
#'
#' @return the vector of consensus clustering result
#'
#' @details
#' Perform majority voting as a consensus mechanism.
#'
#' @examples
#' X = gaussian_clusters()$X
#' x1 = kmeans(X, 5)$cluster
#' x2 = kmeans(X, 5)$cluster
#' x3 = kmeans(X, 5)$cluster
#' clusters = majority_voting(cbind(x1,x2,x3))
#'
majority_voting = function(X){

  Ncusters = ncol(X)

  Clusters = X[,1]
  for (i in 2:Ncusters){
    clusters = cluster_relabel(X[,1], X[,i])
    Clusters = cbind(Clusters, clusters[,3])
  }

  clusters = apply(Clusters, 1, function(x){return(names(which.max(table(x))))})
  return(clusters)
}

