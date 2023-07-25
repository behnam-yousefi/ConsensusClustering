## Consensus clustering functions_...
## version that works for multiview data in cases where not all the data is available for all the datasets.
## Author: Behnam Yousefi
library(igraph)
library(cluster)

## Basic Functions -------------------------------------------------------------
Logit = function(x)
  return(log(1*x/(1-x)))

# Build connectivity matrix for resampled clustering
ConnMat = function (Clusters){
  Nsample = length(Clusters)
  M = matrix(0,Nsample,Nsample)
  for (i in 1:Nsample)
    for (j in 1:Nsample)
      if (Clusters[i]*Clusters[j]>0)
        M[i,j] = ifelse(Clusters[i]==Clusters[j],1,0)
  return(M)
}

#
IndcMat = function (Clusters){
  Nsample = length(Clusters)
  I = matrix(0,Nsample,Nsample)
  for (i in 1:Nsample)
    for (j in 1:Nsample)
      I[i,j] = ifelse(Clusters[i]*Clusters[j]>0,1,0)
  return(I)
}

# Convert adjacency function to the affinity matrix (Luxburg et al.)
AdjConv = function(Adj, alpha=1){
  # Adj must be within [-1, 1]
  Aff = exp(-(1-Adj)^2/(2*alpha^2))
  return(Aff)
}

# Clustering methods
HirClustFromAdjMat = function (Adj, k = 2, alpha = 1, adj.conv = TRUE, method = "ward.D"){
  Adj[is.na(Adj)] = 0

  Aff = if(adj.conv)
    Aff = AdjConv(Adj, alpha)
  else
    Aff = Adj

  dists = as.dist(1-Aff)
  Tree = hclust(dists, method=method)
  clusters = cutree(Tree, k=k)

  return(clusters)
}

SpectClustFromAdjMat = function (Adj, k = 2, max.eig = 10, alpha = 1, adj.conv = TRUE, do.plot = FALSE){
  Adj[is.na(Adj)] = 0

  Aff = if(adj.conv)
    Aff = AdjConv(Adj, alpha)
  else
    Aff = Adj

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

# Clustering methods
PamClustFromAdjMat = function (Adj, k = 2, alpha = 1, adj.conv = TRUE){
  Adj[is.na(Adj)] = 0

  Aff = if(adj.conv)
    Aff = AdjConv(Adj, alpha)
  else
    Aff = Adj

  dists = as.dist(1-Aff)
  clusters = pam(as.matrix(dists), k = k, diss = TRUE, cluster.only=TRUE)

  return(clusters)
}

## Co-clustering matrix -----------------------------------------------------

coCluster_matrix = function(X){
  # Calculate the  Co-cluster matrix for a set of given set of clustering results.
  # Clusters: matrix of Nsamples x Nclusterings
  # Output: The normalized matrix of Co-cluster frequency of any pairs of samples (Nsamples x Nsamples)
  # Note that zeros are are considered as unclustered samples

  Nsample = nrow(X)
  Nmethod = ncol(X)

  M = matrix(0,Nsample,Nsample)
  I = M
  pb = txtProgressBar(min = 0, max = Nmethod, style = 3)

  for (cl in 1:Nmethod){
    M = M + ConnMat(X[,cl])
    I = I + IndcMat(X[,cl])
    setTxtProgressBar(pb, cl)
  }
  coClusterMatrix = M/I
  close(pb)

  rownames(coClusterMatrix) = rownames(X)
  colnames(coClusterMatrix) = rownames(X)
  return(coClusterMatrix)
}

## Multiple K-means -----------------------------------------------------
# 1. Same data
multi_kmeans = function(X, rep = 10, range.k = c(2,5), method = "random"){
  # X: Sample x feature matrix
  # K is selected randomly from a discrete uniform distribution between range.k[1] and range.k[2]
  # method = "random", "silhouette"

  assertthat::assert_that(rep > 0)
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] >= range.k[1])

  Kmin = range.k[1]
  Kmax = range.k[2]
  distX = dist(X)

  Clusters = matrix(0, nrow(X), rep)
  for (i in 1:rep){

    if (method == "silhouette"){
      Sil = rep(0, Kmax)
      for (k in Kmin:Kmax){
        clusters = kmeans(X, k)$cluster
        sil = silhouette(clusters, distX)
        Sil[k] = mean(sil[,"sil_width"])
      }
      Kopt = which.max(Sil)

    }
    else if (method == "random"){
      Kopt = sample(Kmin:Kmax, 1)
    }
    else
      error("err")

    Clusters[,i] = kmeans(X, Kopt)$cluster
    print(Kopt)
  }

  return(Clusters)
}

# 2. Different data
multiview_kmeans = function(X, rep = 10, range.k = c(2,5), method = "random"){
  # X: List of Sample x feature matrices
  # K is selected based on mean Silhouette index for each element of X

  assertthat::assert_that(is.list(X))
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] >= range.k[1])

  Kmin = range.k[1]
  Kmax = range.k[2]
  Nview = length(X)

  Clusters = c()
  for (i in 1:Nview){
    X_i = X[[i]]
    cl = multi_kmeans (X_i, rep = rep, range.k = range.k, method = method)
    Clusters = cbind(Clusters, cl)
  }

  return(Clusters)
}

## Multiple PAM -----------------------------------------------------
# 1. Same data
multi_pam = function(X, rep = 10, range.k = c(2,5), is.distance = FALSE, method = "random"){
  # X: Sample x feature matrix
  # K is selected randomly from a discrete uniform distribution between range.k[1] and range.k[2]
  # method = "random", "silhouette"

  assertthat::assert_that(rep > 0)
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] >= range.k[1])

  Kmin = range.k[1]
  Kmax = range.k[2]

  Clusters = matrix(0, nrow(X), rep)
  for (i in 1:rep){

    if (method == "silhouette"){
      Sil = rep(0, Kmax)
      for (k in Kmin:Kmax){
        pam_result = pam(X, k = k, diss = is.distance)
        clusters = pam_result$clustering
        Sil[k] = pam_result$silinfo$avg.width
      }
      Kopt = which.max(Sil)

    }
    else if (method == "random"){
      Kopt = sample(Kmin:Kmax, 1)
    }
    else
      error("err")

    Clusters[,i] = pam(X, k = Kopt, diss = is.distance, cluster.only = TRUE)
    print(Kopt)
  }

  return(Clusters)
}

# 2. Different data
multiview_pam = function(X, rep = 10, range.k = c(2,5), is.distance = FALSE, method = "random", sample.set = NA){
  # X: List of Sample x feature matrices
  # K is selected based on mean Silhouette index for each element of X
  ## sample.set: a set of samples the clustering is being applied on. can be names or indices
  ## if sample.set is NA, we consider all the datasets have the same samples with the same order
  # method = "random", "silhouette"

  assertthat::assert_that(is.list(X))
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] >= range.k[1])

  if (!is.na(sample.set)[1] & is.null(colnames(X[[1]])))
    error("err")

  Kmin = range.k[1]
  Kmax = range.k[2]
  Nview = length(X)

  Clusters = c()
  for (i in 1:Nview){
    X_i = X[[i]]

    if (!is.na(sample.set)[1]){

      IntersectedSamples = intersect(sample.set, colnames(X_i))

      if (is.distance)
        X_i = X_i[IntersectedSamples,IntersectedSamples]
      else
        X_i = X_i[,IntersectedSamples]

    for (r in 1:rep){
      cl = multi_pam (X_i, rep = 1, range.k = range.k, is.distance = is.distance, method = method)
      names(cl) = IntersectedSamples
      clusters = rep(0, length(sample.set))
      names(clusters) = sample.set
      clusters[names(cl)] = cl
      Clusters = cbind(Clusters, clusters)
    }

    }else{
      cl = multi_pam (X_i, rep = rep, range.k = range.k, is.distance = is.distance, method = method)
      Clusters = cbind(Clusters, cl)
    }
  }
  return(Clusters)
}

## consensus matrix calculation for each Nclust -------------------------------------------
# 1. Same data
consensus_matrix = function(X, max.cluster = 5, resample.ratio = 0.7, max.itter = 100, clustering.method = "hclust",
                            no.cores = 1, adj.conv = TRUE){
  ## Monti et al. (2003) consensus clustering algorithm
  ## X is the adjacency matrix a Nsample x Nsample matrix
  ## The output is the consensus matrix for each k
  ## method: "hclust", "spectral", "pam"

  assertthat::assert_that(ncol(X) == nrow(X))
  assertthat::assert_that(max.cluster>2)
  assertthat::assert_that(resample.ratio>0 & resample.ratio<1)
  assertthat::assert_that(max.itter>0)
  assertthat::assert_that(no.cores>0)

  Nsample = nrow(X)
  CM = list()      # List of consensus matrices

  print("Algorithm Starts...!")
  for (nClust in 2:max.cluster) {                          # Loop for K

    print(paste0("Number of clusters: ", nClust))

    ## Connectivity matrix
    M = matrix(0,Nsample,Nsample)
    rownames(M) = rownames(X)
    colnames(M) = rownames(X)
    ## Indicator matrix
    I = M

    pb = txtProgressBar(min = 0, max = max.itter, style = 3)
    for (i in 1 : max.itter){
      setTxtProgressBar(pb, i)

      RandInd = sample(Nsample, floor(resample.ratio*Nsample), replace = FALSE)
      X_i = X[RandInd,RandInd]

      ## Do clustering
      if (clustering.method == "hclust")
        clusters = HirClustFromAdjMat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "spectral")
        clusters = SpectClustFromAdjMat(X_i, k = nClust, max.eig = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "pam")
        clusters = PamClustFromAdjMat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else
        error("err")
      Clusters = rep(0,Nsample)
      names(Clusters) = rownames(X)
      Clusters[RandInd] = clusters

      ## Connectivity matrix
      Mi = ConnMat(Clusters)
      M = M + Mi

      Ii = IndcMat(Clusters)
      I = I + Ii
    }
    close(pb)

    CM[[nClust]] = M/I
    rownames(CM[[nClust]]) = rownames(X)
    colnames(CM[[nClust]]) = rownames(X)
  }

  return(CM)
}

# 1. Same data
multiview_consensus_matrix = function(X, max.cluster = 5, sample.set = NA, clustering.method = "hclust",
                                      no.cores = 1, adj.conv = TRUE){
  ## Monti et al. (2003) consensus clustering algorithm
  ## X is a list of different data (or view) each a Nsample x Nsample matrix
  ## The output is the consensus matrix for each k
  ## method: "hclust", "spectral", "pam"
  ## sample.set: a set of samples the clustering is being applied on. can be names or indices
  ## if sample.set is NA, we consider all the datasets have the same samples with the same order

  assertthat::assert_that(is.list(X))
  assertthat::assert_that(max.cluster>=2)
  assertthat::assert_that(no.cores>0)

  if (is.na(sample.set)[1])
    sample.set = 1:ncol(X[[1]])

  N_dataset = length(X)
  CM = list()      # List of consensus matrices

  print("Algorithm Starts...!")
  for (nClust in 2:max.cluster) {                          # Loop for K

    print(paste0("Number of clusters: ", nClust))

    ## Connectivity matrix
    M = matrix(0,length(sample.set),length(sample.set))
    rownames(M) = sample.set
    colnames(M) = sample.set
    ## Indicator matrix
    I = M

    pb = txtProgressBar(min = 0, max = N_dataset, style = 3)
    for (i in 1 : N_dataset){
      setTxtProgressBar(pb, i)

      X_i = X[[i]]
      IntersectedSamples = intersect(sample.set, rownames(X_i))
      X_i = X_i[IntersectedSamples,IntersectedSamples]

      ## Do clustering
      if (clustering.method == "hclust")
        clusters = HirClustFromAdjMat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "spectral")
        clusters = SpectClustFromAdjMat(X_i, k = nClust, max.eig = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "pam")
        clusters = PamClustFromAdjMat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else
        error("err")
      Clusters = rep(0,length(sample.set))
      names(Clusters) = sample.set
      Clusters[names(clusters)] = clusters

      ## Conectivity matrix
      Mi = ConnMat(Clusters)
      M = M + Mi

      Ii = IndcMat(Clusters)
      I = I + Ii
    }
    close(pb)

    CM[[nClust]] = M/I
    rownames(CM[[nClust]]) = rownames(X)
    colnames(CM[[nClust]]) = rownames(X)
  }

  return(CM)
}

## Cluster count function based on consensus matrices  -------------------------------------------

CC_cluster_count = function(CM, plot.cdf = TRUE, plot.logit = FALSE){

  Nsample = ncol(CM[[2]])
  K = 2:(length(CM))

  par(new=FALSE)
  A = rep(0,length(K))                       # Area under the CDF curve
  RobScore = rep(0,length(K))
  LogitScore = rep(0,length(K))
  PAC = rep(0,length(K))
  CMavg = rep(0,length(K))

  Nbin = 100

  for (k in K){

    M0 = CM[[k]]

    CMavg[k] = mean(M0)

    diag(M0) = NA

    # CDF
    ConsDistr = rep(0,Nbin-1)
    for (i in 1:(Nbin-1))
      ConsDistr[i] = sum(M0<(i/Nbin), na.rm = TRUE)/2
    ConsDistr = ConsDistr/(Nsample*(Nsample-1)/2)
    A[k] = sum(ConsDistr)/Nbin

    S_Window_L = sum(ConsDistr*dnorm(0:(Nbin-2),10,5))
    S_Window_H = sum(ConsDistr*dnorm(0:(Nbin-2),90,5))
    RobScore[k] = (S_Window_L - 0) * (1 - S_Window_H) / (S_Window_H - S_Window_L)^2

    # logit score
    df = data.frame(x = seq(0.01, .99, length.out = length(ConsDistr)), y = ConsDistr)
    linear_model = lm(y ~ Logit(x), data=df)
    b_hat = linear_model[["coefficients"]][1]
    a_hat = linear_model[["coefficients"]][2]
    LogitScore[k] = a_hat

    S_Window_L = ConsDistr[10]
    S_Window_H = ConsDistr[90]
    PAC[k] = S_Window_H - S_Window_L

    if (plot.cdf){
      # ConsDistr = c(0, ConsDistr, 1)
      plt = ConsDistr
      index_value = seq(0, 1, length.out = length(ConsDistr))
      plot(index_value,ConsDistr, type='l', col=k, ylim = c(0,1), lwd = 3,
           ylab = "CDF", xlab = "consensus index value")
      text(k/10, plt[10*k], labels = k)
      par(new=TRUE)
    }

    if (plot.logit){
      x = seq(.001,.999,.001)
      print(paste0("k = ", k, ", a_hat = ", a_hat))
      y = a_hat * Logit(x) + b_hat
      #x = c(0, x, 1)
      #y = c(0, y, 1)
      plot(x, y, type='l', col=k, ylim = c(0,1), lwd = 3,
           ylab = "simulated CDF", xlab = "consensus index value")
      text(k/10, ConsDistr[10*k], labels = k)
      par(new=TRUE)
    }

  }
  
  par(new=FALSE)

  deltaA = A
  for (k in 3:max(K))
    deltaA[k] = (A[k]-A[k-1])/A[k-1]

  Kopt_LogitScore = which.min(LogitScore[-1]) + 1
  if (Kopt_LogitScore == max(K)){
    # look at the knee point
    LogitScore_Knee = rep(0, max(K))
    for (k in 3:(max(K)-1)){
      m1 = LogitScore[k] - LogitScore[k-1]
      m2 = LogitScore[k+1] - LogitScore[k]
      LogitScore_Knee[k] = m2 - m1
    }
    # plot(LogitScore_Knee)
    Kopt_LogitScore = which.max(LogitScore_Knee)
  }

  Result = list()
  Result[["RobScore"]] = RobScore
  Result[["LogitScore"]] = LogitScore
  Result[["PAC"]] = PAC
  Result[["deltaA"]] = deltaA
  Result[["CMavg"]] = CMavg
  Result[["Kopt_RobScore"]] = which.max(RobScore)
  Result[["Kopt_LogitScore"]] = Kopt_LogitScore
  Result[["Kopt_PAC"]] = which.min(PAC[-1]) + 1
  Result[["Kopt_deltaA"]] = which.max(deltaA)
  Result[["Kopt_CMavg"]] = which.max(CMavg)
  return(Result)
}

## Null distribution of robust score generation by permuting the Adj -------------------------------------------

robustness_null_dist = function(Adj, rep = 10, max.cluster = 5, resample.ratio = 0.7, max.itter = 100, clustering.method = "hclust",
                                 is.similarity = TRUE, no.cores = 1, adj.conv = TRUE){

  print("Calculation of the null distribution of robust score ...")

  ScoreList = list()
  ScoreList[["RobScore"]] = matrix(0, rep, max.cluster)
  ScoreList[["LogitScore"]] = matrix(0, rep, max.cluster)
  ScoreList[["PAC"]] = matrix(0, rep, max.cluster)
  ScoreList[["deltaA"]] = matrix(0, rep, max.cluster)
  ScoreList[["CMavg"]] = matrix(0, rep, max.cluster)

  N = ncol(Adj)

  for (i in 1:rep){
    Adj_permuted = matrix(sample(Adj, N*N, replace = TRUE), N,N)
    # Adj_permuted = SampelGraph_Degree(Adj)
    CM = consensus_matrix(Adj_permuted, max.cluster = max.cluster, resample.ratio = resample.ratio,
                          max.itter = max.itter, clustering.method = clustering.method, adj.conv = adj.conv)
    Scores = CC_cluster_count(CM)
    ScoreList[["RobScore"]][i,] = Scores[["RobScore"]]
    ScoreList[["LogitScore"]][i,] = Scores[["LogitScore"]]
    ScoreList[["PAC"]][i,] = Scores[["PAC"]]
    ScoreList[["deltaA"]][i,] = Scores[["deltaA"]]
    ScoreList[["CMavg"]][i,] = Scores[["CMavg"]]
  }

  return(ScoreList)
}

## Calculate the null-model corrected score -------------------------------------------

corrected_score = function(scores, null_model_scores, method = "log.ratio"){

  MeanNull = apply(null_model_scores, 2, mean)
  corrected_score =  abs(log10(scores + 1e-10) - log10(MeanNull + 1e-10))

  return(corrected_score)

}

## Calculate the null-model p-values -------------------------------------------

cluster_Pvalues = function(scores, null_model_scores){

  n_total = nrow(null_model_scores)

  p_val = rep(1, length(scores))
  for (i in 1:length(scores)){
    n_less = sum(scores[i] > null_model_scores[,i])
    p_val[i] = (n_less + 1) / (n_total + 1)
  }

  return(p_val)
}

## Graph sampling with constant node degree
SampelGraph_Degree = function(A){

  EdgeWeight = matrix(A, length(A),1)
  EdgeWeight = EdgeWeight[EdgeWeight>0]

  Nnode = ncol(A)
  Degree = apply(A, 1, sum)
  As = A
  for (i in 1:Nnode-1){
    Weights = sample(EdgeWeight,Nnode-i)
    Weights = Degree[i] * Weights / (sum(Weights) + sum(As[i,1:i]))
    As[i,(i+1):Nnode] = Weights
    As[,i] = As[i,]
  }
  return(As)
}
