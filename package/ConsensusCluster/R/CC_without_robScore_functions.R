# Generation mechanisms without robustness score
# Functions accept a specific range of hyperparameters (specifically number of clusters)
# and perform several clusterings. The output, that is a matrix of clusterings, can then
# be fed to a consensus function to obtain the final clustering.

# ---------------------- Multi-Method Generation -----------------------

#' Multiple K-means generation
#'
#' @param X input data Nsample x Nfeatures
#' @param rep number of repeats
#' @param range.k vector of minimum and maximum values for k \code{c(min, max)}
#' @param method method for the choice of k at each repeat \code{c("random", "silhouette")}
#'
#' @return matrix of clusterings Nsample x Nrepeat
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then k-means clustering is applied and result is returned.
#'
#' @examples
#' X = gaussian_clusters()$X
#' Clusters = multi_kmeans_gen(X)
#'
multi_kmeans_gen = function(X, rep = 10, range.k = c(2,5), method = "random"){

  assertthat::assert_that(rep > 0)
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] >= range.k[1])
  assertthat::assert_that(method %in% c("silhouette", "random"))

  Kmin = range.k[1]
  Kmax = range.k[2]
  distX = dist(X)

  Clusters = matrix(0, nrow(X), rep)
  for (i in 1:rep){

    if (method == "silhouette"){
      Sil = rep(0, Kmax)
      for (k in Kmin:Kmax){
        clusters = stats::kmeans(X, k)$cluster
        sil = cluster::silhouette(clusters, distX)
        Sil[k] = mean(sil[,"sil_width"])
      }
      Kopt = which.max(Sil)

    }
    else if (method == "random"){
      Kopt = sample(Kmin:Kmax, 1)
    }
    else
      error("err")

    Clusters[,i] = stats::kmeans(X, Kopt)$cluster
    # print(Kopt)
  }

  return(Clusters)
}


#' Multiple PAM (K-medoids) generation
#'
#' @param X input data Nsample x Nfeatures or distance matrix.
#' @param rep number of repeats
#' @param range.k vector of minimum and maximum values for k \code{c(min, max)}
#' @param is.distance binary balue indicating if the input \code{X} is distance
#' @param method method for the choice of k at each repeat \code{c("random", "silhouette")}
#'
#' @return matrix of clusterings Nsample x Nrepeat
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then PAM clustering is applied and result is returned.
#'
#' @examples
#' X = gaussian_clusters()$X
#' Clusters = multi_pam_gen(X)
#'
multi_pam_gen = function(X, rep = 10, range.k = c(2,5), is.distance = FALSE, method = "random"){

  assertthat::assert_that(rep > 0)
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] >= range.k[1])
  assertthat::assert_that(method %in% c("silhouette", "random"))

  Kmin = range.k[1]
  Kmax = range.k[2]

  Clusters = matrix(0, nrow(X), rep)
  for (i in 1:rep){

    if (method == "silhouette"){
      Sil = rep(0, Kmax)
      for (k in Kmin:Kmax){
        pam_result = cluster::pam(X, k = k, diss = is.distance)
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

    Clusters[,i] = cluster::pam(X, k = Kopt, diss = is.distance, cluster.only = TRUE)
    # print(Kopt)
  }

  return(Clusters)
}

#' Multiple cluster generation
#'
#' @param X input data Nsample x Nfeatures or a distance matrix
#' @param func custom function that accepts \code{X} and a parameter that return a vector of clusterings.
#' \code{cluster_func <- function(X, param)}
#' @param rep number of repeats
#' @param param vector of parameters
#' @param method method for the choice of k at each repeat \code{c("random", "silhouette")}
#'
#' @return matrix of clusterings Nsample x Nrepeat
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then clustering is applied and result is returned.
#'
#' @examples
#' X = gaussian_clusters()$X
#' cluster_func = function(X, k){return(stats::kmeans(X, k)$cluster)}
#' Clusters = multi_cluster_gen(X, cluster_func, param = c(2,3))
#'
#'
multi_cluster_gen = function(X, func, rep = 10, param, method = "random"){

  assertthat::assert_that(rep > 0)
  assertthat::assert_that(method %in% c("silhouette", "random"))

  distX = dist(X)

  Clusters = matrix(0, nrow(X), rep)
  for (i in 1:rep){

    if (method == "silhouette"){

      Sil = rep(0, length(param))
      for (k in param){
        clusters = func(X, k)
        sil = cluster::silhouette(clusters, distX)
        Sil[k] = mean(sil[,"sil_width"])
      }
      Kopt = which.max(Sil)

    }
    else if (method == "random"){
      Kopt = sample(param, 1)
    }
    else
      error("err")

    Clusters[,i] = func(X, Kopt)
    # print(Kopt)
  }

  return(Clusters)
}


# ---------------------- Multi-Data CC -----------------------

#' Multiview K-means generation
#'
#' @param X List of input data matrices of Sample x feature. The length of \code{X} is equal to Nviews
#' @param rep number of repeats
#' @param range.k vector of minimum and maximum values for k \code{c(min, max)}
#' @param method method for the choice of k at each repeat \code{c("random", "silhouette")}
#'
#' @return matrix of clusterings Nsample x (Nrepeat x Nviews)
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then k-means clustering is applied and result is returned.
#'
#' @examples
#' multiview_kmeans_gen(X)
#'
#'
multiview_kmeans_gen = function(X, rep = 10, range.k = c(2,5), method = "random"){

  assertthat::assert_that(is.list(X))
  assertthat::assert_that(range.k[1] > 1)
  assertthat::assert_that(range.k[2] >= range.k[1])

  Kmin = range.k[1]
  Kmax = range.k[2]
  Nview = length(X)

  Clusters = c()
  for (i in 1:Nview){
    X_i = X[[i]]
    cl = multi_kmeans_gen(X_i, rep = rep, range.k = range.k, method = method)
    Clusters = cbind(Clusters, cl)
  }

  return(Clusters)
}

#' Multiview PAM (K-medoids) generation
#'
#' @param X List of input data matrices of Sample x feature or distance matrices.
#' The length of \code{X} is equal to Nviews
#' @param rep number of repeats
#' @param range.k vector of minimum and maximum values for k \code{c(min, max)}
#' @param is.distance binary balue indicating if the input \code{X} is distance
#' @param method method for the choice of k at each repeat \code{c("random", "silhouette")}
#' @param sample.set  vector of samples the clustering is being applied on. can be names or indices.
#' if \code{sample.set} is \code{NA}, it considers all the datasets have the same samples with the same order
#'
#' @return matrix of clusterings Nsample x (Nrepeat x Nviews)
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then PAM clustering is applied and result is returned.
#'
#' @examples
#' multiview_kmeans_gen(X)
#'
#'
multiview_pam_gen = function(X, rep = 10, range.k = c(2,5), is.distance = FALSE, method = "random", sample.set = NA){

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
        cl = multi_pam_gen(X_i, rep = 1, range.k = range.k, is.distance = is.distance, method = method)
        names(cl) = IntersectedSamples
        clusters = rep(0, length(sample.set))
        names(clusters) = sample.set
        clusters[names(cl)] = cl
        Clusters = cbind(Clusters, clusters)
      }

    }else{
      cl = multi_pam_gen(X_i, rep = rep, range.k = range.k, is.distance = is.distance, method = method)
      Clusters = cbind(Clusters, cl)
    }
  }
  return(Clusters)
}

#' Multiview cluster generation
#'
#' @param X List of input data matrices of Sample x feature or distance matrices.
#' The length of \code{X} is equal to Nviews
#' @param func custom function that accepts \code{X} and a parameter that return a vector of clusterings.
#' \code{cluster_func <- function(X, param)}
#' @param rep number of repeats
#' @param param vector of parameters
#' @param method method for the choice of k at each repeat \code{c("random", "silhouette")}
#' @param is.distance binary balue indicating if the input \code{X[i]} is distance
#' @param sample.set  vector of samples the clustering is being applied on. can be names or indices.
#' if \code{sample.set} is \code{NA}, it considers all the datasets have the same samples with the same order
#'
#' @return matrix of clusterings Nsample x (Nrepeat x Nviews)
#'
#' @details
#' At each repeat, k is selected randomly or based on the best silhouette width from a discrete uniform distribution between range.k[1] and range.k[2].
#' Then clustering is applied and result is returned.
#'
#' @examples
#' multiview_kmeans_gen(X)
#'
#'
multiview_pam_gen = function(X, func, rep = 10, param, method = "random", is.distance = FALSE, sample.set = NA){

  assertthat::assert_that(is.list(X))

  if (!is.na(sample.set)[1] & is.null(colnames(X[[1]])))
    error("err")

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
        cl = multi_cluster_gen(X_i, rep = 1, param = param, method = method)
        names(cl) = IntersectedSamples
        clusters = rep(0, length(sample.set))
        names(clusters) = sample.set
        clusters[names(cl)] = cl
        Clusters = cbind(Clusters, clusters)
      }

    }else{
      cl = multi_cluster_gen(X_i, rep = 1, param = param, method = method)
      Clusters = cbind(Clusters, cl)
    }
  }
  return(Clusters)
}
