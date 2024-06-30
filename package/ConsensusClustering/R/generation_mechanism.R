# Generation mechanism functions

#' Generation mechanism for data perturbation consensus clustering
#'
#' @param X input data Nsample x Nfeatures
#' @param cluster.method base clustering method: \code{c("hclust", "spectral", "pam", "custom")}
#' @param k number of clusters
#' @param resample.ratio the data ratio to use at each itteration.
#' @param rep maximum number of itterations at each \code{max.cluster}
#' @param distance.method method for distance calculation:
#' \code{"euclidian"}, \code{"cosine"}, \code{"maximum"}, \code{"manhattan"},
#' \code{"canberra"}, \code{"binary"}, \code{"minkowski"}.
#' @param adj.conv binary value to apply soft thresholding (default=\code{TRUE})
#' @param func user-definrd function required if \code{cluster.method = "custom"}.
#' The function needs two inputs of X and k
#'
#' @return matrix of clusterings Nsample x Nrepeat
#'
#' @details
#' Performs clustering on the purturbed samples set
#' Monti et al. (2003) consensus clustering algorithm
#'
#' @examples
#' X = gaussian_clusters()$X
#' Clusters = generate_data_prtrb(X)
#'
generate_data_prtrb = function(X, cluster.method = "pam", k = 3, resample.ratio = 0.7,
                               rep = 10, distance.method = "euclidian", adj.conv = TRUE, func){

  assertthat::assert_that(k>2)
  assertthat::assert_that(resample.ratio>0 & resample.ratio<1)
  assertthat::assert_that(rep > 0)

  Nsample = nrow(X)

  # Convert X into distance matrix
  X = adj_mat(X, method = distance.method)

  # Perform the repeat loop
  Clusters = matrix(0, Nsample, rep)
  for (i in 1:rep){
    ## Sample randomization
    RandInd = sample(Nsample, floor(resample.ratio*Nsample), replace = FALSE)
    X_i = X[RandInd,RandInd]
    # sd.noise = .01
    # if (sd.noise > 0)
    #   X_i = X_i + matrix(rnorm(nrow(X_i)*ncol(X_i), 0, sd.noise), nrow(X_i), ncol(X_i))

    ## Do clustering
    if (cluster.method == "hclust")
      clusters = hir_clust_from_adj_mat(X_i, k = k, alpha = 1, adj.conv = adj.conv)
    else if (cluster.method == "spectral")
      clusters = spect_clust_from_adj_mat(X_i, k = k, max.eig = k, alpha = 1, adj.conv = adj.conv)
    else if (cluster.method == "pam")
      clusters = pam_clust_from_adj_mat(X_i, k = k, alpha = 1, adj.conv = adj.conv)
    else if (cluster.method == "custom")
      clusters = func(X_i, k)
    else
      stop("cluster.method not implemented")

    ## Add zeros for absent samples
    clusters_with_0 = rep(0,Nsample)
    names(clusters_with_0) = rownames(X)
    clusters_with_0[RandInd] = clusters
    Clusters[,i] = clusters_with_0
  }
  return(Clusters)
}

#' Multiple method generation
#'
#' @param X input data Nsample x Nfeatures
#' @param cluster.method base clustering method: \code{c("kmeans", "pam", "custom")}
#' @param range.k vector of minimum and maximum values for k \code{c(min, max)}
#' @param sample.k.method method for the choice of k at each repeat \code{c("random", "silhouette")}
#' @param rep number of repeats
#' @param distance.method method for distance calculation:
#' \code{"euclidian"}, \code{"maximum"}, \code{"manhattan"},
#' \code{"canberra"}, \code{"binary"}, \code{"minkowski"}.
#' @param func user-definrd function required if \code{cluster.method = "custom"}.
#' The function needs two inputs of X and k.
#'
#' @return matrix of clusterings Nsample x Nrepeat
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then clustering is applied and result is returned.
#'
#' @examples
#' X = gaussian_clusters()$X
#' Clusters = generate_method_prtrb(X)
#'
generate_method_prtrb = function(X, cluster.method = "pam", range.k = c(2, 5), sample.k.method = "random",
                                 rep = 10, distance.method = "euclidian", func){

  assertthat::assert_that(rep > 0)
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] > range.k[1])
  assertthat::assert_that(sample.k.method %in% c("silhouette", "random"))

  Kmin = range.k[1]
  Kmax = range.k[2]
  Nsample = nrow(X)

  # Perform the repeat loop
  Clusters = matrix(0, Nsample, rep)
  for (i in 1:rep){

    ## Get the vakue of k
    if (sample.k.method == "silhouette"){
      distX = stats::dist(X, distance.method)
      Sil = rep(0, Kmax)
      for (k in Kmin:Kmax){
        if (cluster.method == "kmeans"){
          clusters = stats::kmeans(X, k)$cluster
          sil = cluster::silhouette(clusters, distX)
          Sil[k] = mean(sil[,"sil_width"])
        }else if (cluster.method == "pam"){
          pam_result = cluster::pam(X, k = k, diss = FALSE)
          clusters = pam_result$clustering
          Sil[k] = pam_result$silinfo$avg.width
        }else if (cluster.method == "custom"){
          clusters = func(X, k)
          sil = cluster::silhouette(clusters, distX)
          Sil[k] = mean(sil[,"sil_width"])
        }else{
          stop("cluster.method not implemented")
        }
      Kopt = which.max(Sil)
      }
    }else if (sample.k.method == "random"){
      Kopt = sample(Kmin:Kmax, 1)
    }else{
      stop("sample.k.method not implemented")
    }

    ##Do cluster
    if (cluster.method == "kmeans")
      Clusters[,i] = stats::kmeans(X, Kopt)$cluster
    else if (cluster.method == "pam")
      Clusters[,i] = cluster::pam(X, k = Kopt, diss = FALSE, cluster.only = TRUE)
    else if (cluster.method == "custom")
      Clusters[,i] = func(X, Kopt)
    else
      stop("cluster.method not implemented")
  }
  return(Clusters)
}

#' Multiview generation
#'
#' @param X list of input data matrices of Sample x feature or distance matrices.
#' The length of \code{X} is equal to Nviews
#' @param cluster.method base clustering method: \code{c("kmeans", "pam", "custom")}
#' @param range.k vector of minimum and maximum values for k \code{c(min, max)}
#' @param sample.k.method method for the choice of k at each repeat \code{c("random", "silhouette")}
#' @param rep number of repeats
#' @param distance.method method for distance calculation:
#' \code{"euclidian"}, \code{"maximum"}, \code{"manhattan"},
#' \code{"canberra"}, \code{"binary"}, \code{"minkowski"}.
#' @param sample.set  vector of samples the clustering is being applied on. can be names or indices.
#' If \code{sample.set} is \code{NA}, it considers all the datasets have the same samples with the same order
#' @param func user-definrd function required if \code{cluster.method = "custom"}.
#' The function needs two inputs of X and k.
#'
#' @return matrix of clusterings Nsample x Nrepeat
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then clustering is applied and result is returned.
#'
#' @examples
#' data = multiview_clusters (n = c(40,40,40), hidden.dim = 2, observed.dim = c(2,2,2),
#' sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))
#' X_observation = data[["observation"]]
#' Clusters = multiview_pam_gen(X_observation)
#'
generate_multiview = function(X, cluster.method = "pam", range.k = c(2, 5), sample.k.method = "random",
                              rep = 10, distance.method = "euclidian", sample.set = NA, func){

  assertthat::assert_that(is.list(X))
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] > range.k[1])
  assertthat::assert_that(sample.k.method %in% c("silhouette", "random"))

  if (!is.na(sample.set)[1] & is.null(colnames(X[[1]])))
    stop("err")

  Kmin = range.k[1]
  Kmax = range.k[2]
  Nview = length(X)

  Clusters = c()
  for (i in 1:Nview){
    X_i = X[[i]]

    if (!is.na(sample.set)[1]){

      IntersectedSamples = intersect(sample.set, colnames(X_i))

      # if (is.distance)
      #  X_i = X_i[IntersectedSamples,IntersectedSamples]
      # else
      X_i = X_i[,IntersectedSamples]

      for (r in 1:rep){
        cl = generate_method_prtrb(X_i, cluster.method = cluster.method, range.k = range.k, sample.k.method = sample.k.method,
                                  rep = 1, distance.method = distance.method, func = func)
        names(cl) = IntersectedSamples
        clusters = rep(0, length(sample.set))
        names(clusters) = sample.set
        clusters[names(cl)] = cl
        Clusters = cbind(Clusters, clusters)
      }

    }else{
      cl = generate_method_prtrb(X_i, cluster.method = cluster.method, range.k = range.k, sample.k.method = sample.k.method,
                                 rep = rep, distance.method = distance.method, func = func)
      Clusters = cbind(Clusters, cl)
    }
  }
  return(Clusters)
}
