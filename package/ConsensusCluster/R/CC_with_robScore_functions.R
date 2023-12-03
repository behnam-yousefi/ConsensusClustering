# Generation mechanisms with robustness score
# Functions accept a range of hyperparameters (specifically number of clusters)
# and choose the most optimum based on a robustness (stability) score.
# The output, that is a consensus matrix for the best param, can then
# be fed to a clustering function to obtain the final clustering.

#' Calculate consensus matrix for data perturbation consensus clustering
#'
#' @param X adjacency matrix a Nsample x Nsample
#' @param max.cluster maximum number of clusters
#' @param resample.ratio the data ratio to use at each itteration.
#' @param max.itter maximum number of itterations at each \code{max.cluster}
#' @param clustering.method base clustering method: \code{c("hclust", "spectral", "pam")}
#' @param adj.conv binary value to apply soft thresholding (default=\code{TRUE})
#' @param verbos binary value for verbosity (default=\code{TRUE})
#'
#' @return list of consensus matrices for each k
#'
#' @details
#' performs data perturbation consensus clustering and obtain consensus matrix
#' Monti et al. (2003) consensus clustering algorithm
#'
#' @examples
#' X = gaussian_clusters()$X
#' Adj = adj_mat(X, method = "euclidian")
#' CM = consensus_matrix(Adj, max.cluster=3, max.itter=10, verbos = FALSE)
#'
consensus_matrix = function(X, max.cluster = 5, resample.ratio = 0.7, max.itter = 100, clustering.method = "hclust",
                            adj.conv = TRUE, verbos = TRUE){

  assertthat::assert_that(ncol(X) == nrow(X))
  assertthat::assert_that(max.cluster>2)
  assertthat::assert_that(resample.ratio>0 & resample.ratio<1)
  assertthat::assert_that(max.itter>0)

  Nsample = nrow(X)
  CM = list()      # List of consensus matrices

  if (verbos)
    print("Algorithm Starts...!")
  for (nClust in 2:max.cluster) {                          # Loop for K

    if (verbos)
      print(paste0("Number of clusters: ", nClust))

    ## Connectivity matrix
    M = matrix(0,Nsample,Nsample)
    rownames(M) = rownames(X)
    colnames(M) = rownames(X)
    ## Indicator matrix
    I = M

    if (verbos)
      pb = utils::txtProgressBar(min = 0, max = max.itter, style = 3)
    for (i in 1 : max.itter){

      if (verbos)
        utils::setTxtProgressBar(pb, i)

      RandInd = sample(Nsample, floor(resample.ratio*Nsample), replace = FALSE)
      X_i = X[RandInd,RandInd]

      ## Do clustering
      if (clustering.method == "hclust")
        clusters = hir_clust_from_adj_mat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "spectral")
        clusters = spect_clust_from_adj_mat(X_i, k = nClust, max.eig = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "pam")
        clusters = pam_clust_from_adj_mat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else
        stop("err")
      Clusters = rep(0,Nsample)
      names(Clusters) = rownames(X)
      Clusters[RandInd] = clusters

      ## Connectivity matrix
      Mi = connectivity_matrix(Clusters)
      M = M + Mi

      Ii = indicator_matrix(Clusters)
      I = I + Ii
    }

    if (verbos)
      close(pb)

    CM[[nClust]] = M/I
    rownames(CM[[nClust]]) = rownames(X)
    colnames(CM[[nClust]]) = rownames(X)
  }

  return(CM)
}

#' Calculate consensus matrix for multi-data consensus clustering
#'
#' @param X list of adjacency matrices for different cohorts (or views).
#' @param max.cluster maximum number of clusters
#' @param sample.set vector of samples the clustering is being applied on. \code{sample.set} can be names or indices.
#' if \code{sample.set} is \code{NA}, it considers that all the datasets have the same samples with the same order.
#' @param clustering.method base clustering method: \code{c("hclust", "spectral", "pam")}
#' @param adj.conv binary value to apply soft threshold (default=\code{TRUE})
#' @param verbos binary value for verbosity (default=\code{TRUE})
#'
#' @return description list of consensus matrices for each k
#'
#' @details
#' performs multi-data consensus clustering and obtain consensus matrix
#' Monti et al. (2003) consensus clustering algorithm
#'
#' @examples
#' data = multiview_clusters (n = c(40,40,40), hidden.dim = 2, observed.dim = c(2,2,2),
#' sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))
#' X_observation = data[["observation"]]
#' Adj = list()
#' for (i in 1:length(X_observation))
#'   Adj[[i]] = adj_mat(X_observation[[i]], method = "euclidian")
#' CM = multiview_consensus_matrix(Adj, max.cluster = 4, verbos = FALSE)
#'
multiview_consensus_matrix = function(X, max.cluster = 5, sample.set = NA, clustering.method = "hclust",
                                      adj.conv = TRUE, verbos = TRUE){

  assertthat::assert_that(is.list(X))
  assertthat::assert_that(max.cluster>=2)

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

    if (verbos)
      pb = utils::txtProgressBar(min = 0, max = N_dataset, style = 3)
    for (i in 1 : N_dataset){

      if (verbos)
        utils::setTxtProgressBar(pb, i)

      X_i = X[[i]]
      IntersectedSamples = intersect(sample.set, rownames(X_i))
      X_i = X_i[IntersectedSamples,IntersectedSamples]

      ## Do clustering
      if (clustering.method == "hclust")
        clusters = hir_clust_from_adj_mat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "spectral")
        clusters = spect_clust_from_adj_mat(X_i, k = nClust, max.eig = nClust, alpha = 1, adj.conv = adj.conv)
      else if (clustering.method == "pam")
        clusters = pam_clust_from_adj_mat(X_i, k = nClust, alpha = 1, adj.conv = adj.conv)
      else
        stop("err")
      Clusters = rep(0,length(sample.set))
      names(Clusters) = sample.set
      Clusters[names(clusters)] = clusters

      ## Conectivity matrix
      Mi = connectivity_matrix(Clusters)
      M = M + Mi

      Ii = indicator_matrix(Clusters)
      I = I + Ii
    }

    if (verbos)
      close(pb)

    CM[[nClust]] = M/I
    rownames(CM[[nClust]]) = rownames(X)
    colnames(CM[[nClust]]) = rownames(X)
  }

  return(CM)
}

#' Count the number of clusters based on stability score.
#'
#' @param CM list of consensus matrices each for a specific number of clusters.
#' It can be the output of \code{consensus_matrix()} and \code{multiview_consensus_matrix()} functions.
#' @param plot.cdf binary value to plot the cumulative distribution functions of \code{CM} (default \code{TRUE}).
#' @param plot.logit binary value to plot the logit model of cumulative distribution functions of \code{CM} (default \code{FALSE}).
#'
#' @return results as a list:
#' \code{"LogitScore", "PAC", "deltaA", "CMavg"},
#' \code{"Kopt_LogitScore", "Kopt_PAC", "Kopt_deltaA", "Kopt_CMavg"}
#'
#' @details
#' Count the number of clusters given a list of consensus matrices each for a specific number of clusters.
#' Using different methods: \code{"LogitScore", "PAC", "deltaA", "CMavg"}
#'
#' @examples
#' X = gaussian_clusters()$X
#' Adj = adj_mat(X, method = "euclidian")
#' CM = consensus_matrix(Adj, max.cluster=3, max.itter=10)
#' Result = CC_cluster_count(CM, plot.cdf=FALSE)
#'
CC_cluster_count = function(CM, plot.cdf = TRUE, plot.logit = FALSE){

  Nsample = ncol(CM[[2]])
  K = 2:(length(CM))

  par(new=FALSE)
  A = rep(0,length(K))                       # Area under the CDF curve
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
  Result[["LogitScore"]] = LogitScore
  Result[["PAC"]] = PAC
  Result[["deltaA"]] = deltaA
  Result[["CMavg"]] = CMavg
  Result[["Kopt_LogitScore"]] = Kopt_LogitScore
  Result[["Kopt_PAC"]] = which.min(PAC[-1]) + 1
  Result[["Kopt_deltaA"]] = which.max(deltaA)
  Result[["Kopt_CMavg"]] = which.max(CMavg)
  return(Result)
}
