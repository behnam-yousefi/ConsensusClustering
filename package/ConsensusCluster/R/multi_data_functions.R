# Multiple clustering with different data (multi-view, multi-cohort)

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
