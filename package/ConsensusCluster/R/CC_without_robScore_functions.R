# Multiple clustering with different hyperparameters of a single dataset

## Multiple K-means -----------------------------------------------------

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
