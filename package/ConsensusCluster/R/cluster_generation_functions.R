## Random data generation functions

#' ??
#'
#' @param x x
#' @param n n
#' @value
#' ???
#' @examples
#' divide_interval_int(10,3)
#'
divide_interval_int = function(x,n){
  y = rep(floor(x/n), n-1)
  y[n] = x - sum(y)
  return(y)
}

#' Generate a set of data points from Gaussian distribution
#'
#' @param n number of generated data points
#' @param center data center of desired dimension
#' @param sigma covariance matrix
#' @param label data labels
#'
#' @value
#' Generated data points from Gaussian distribution with given parameters
#'
#' @examples
#' generate_gaussian_data(10, center=c(0,0), sigma=diag(c(1,1)), label=rep(1,10))
#'
generate_gaussian_data = function(n, center=0, sigma=1, label=NA) {
  data = mvtnorm::rmvnorm(n, mean = center, sigma = sigma)
  data = data.frame(data)
  if (!is.na(label))
    data = data %>% dplyr::mutate(class=factor(label))
  return(data)
}

#' Generate clusters of data points from Gaussian distribution with given parameters
#'
#' @param n vector of number of data points in each cluster
#' The length of \code{n} should be equal to the number of clusters.
#' @param center matrix of centers Ncluster x dim
#' @param sigma list of covariance matrices dim X dim.
#' The length of sigma should be equal to the number of clusters.
#'
#' @value
#' matrix of Nsamples x (dim + 1). The last column is cluster labels.
#'
#' @examples
#' center = rbind(c(0,0),
#'                c(1,1))
#' sigma = list(diag(c(1,1)),
#'              diag(2,2))
#' gaussian_clusters_with_param(c(10, 10), center, sigma)
#'
gaussian_clusters_with_param = function(n, center, sigma){

  Ncluster = length(n)
  dim = ncol(center)
  assertthat::assert_that(Ncluster == nrow(center))
  assertthat::assert_that(Ncluster == length(sigma))

  data = c()
  for (i in 1:Ncluster){
    n_i = n[i]
    center_i = center[i,]
    sigma_i = sigma[[i]]
    data_i = generate_gaussian_data(n_i, center_i, sigma_i, i)
    data = rbind(data, data_i)
  }
  return(data)
}

#' Generate clusters of data points from Gaussian distribution with randomly generated parameters
#'
#' @param n vector of number of data points in each cluster
#' The length of \code{n} should be equal to the number of clusters.
#' @param dim number of dimensions
#' @param sd.max maximum standard deviation of clusters
#' @param sd.noise standard deviation of the added noise
#' @param r.range the range (min, max) of distance of cluster centers from the origin
#'
#' @value
#' a list of data points (X) and cluster labels (class)
#'
#' @examples
#' data = gaussian_clusters()
#' X = data$X
#' y = data$class
#'
gaussian_clusters = function(n = c(50,50), dim = 2, sd.max = .1, sd.noise = .01, r.range = c(.1,1)){
  # k: number of clusters
  # n: number of datapoints

  k = length(n)
  assertthat::assert_that(k>=2)

  # Calculate angles by dividing the circule in k equal parts
  # Calculate centers as points within a sphere of radius in [r.min, r.max]
  Angle = seq(0,2*pi,2*pi/k)[1:k]
  Radius = runif(k, r.range[1], r.range[2])

  Centers = matrix(0, k, dim)
  Sigma = list()
  for (i in 1:k){
    theta = Angle[i]
    r = Radius[i]
    for (j in 1:(dim-1))
      Centers[i,j] = r * sin(theta) * cos(theta)^(j-1)
    # old one: Centers[i,dim] = r * sign(pi - theta - .0001)^dim * cos(theta)^(dim-1)
    Centers[i,dim] = r * cos(theta)^(dim-1)

    # Set random covariance matrices
    sigma = matrix(runif(dim*dim, -sd.max/10, sd.max/10), dim)
    sigma[lower.tri(sigma)] = t(sigma)[lower.tri(sigma)]
    diag(sigma) = runif(dim, sd.max/10, sd.max)
    Sigma[[i]] = sigma
  }

  data = gaussian_clusters_with_param(n = n, center = Centers, sigma =  Sigma)
  data[,1:dim] = data[,1:dim] + matrix(rnorm(sum(n)*dim, 0, sd.noise), sum(n))

  data = list(X = data[, 1:dim], class = data[, dim+1])
  return(data)
}

#' Generate clusters of data points from Gaussian-mixture-model distributions with randomly generated parameters
#'
#' @param n vector of number of data points in each cluster
#' The length of \code{n} should be equal to the number of clusters.
#' @param dim number of dimensions
#' @param sd.max maximum standard deviation of clusters
#' @param sd.noise standard deviation of the added noise
#' @param r.range the range (min, max) of distance of cluster centers from the origin
#' @param mixture.range range (min, max) of the number of Gaussian-mixtures.
#' @param mixture.sep scaler indicating the separability between the mixtures.
#'
#' @value
#' a list of data points (X) and cluster labels (class)
#'
#' @examples
#' data = gaussian_mixture_clusters()
#' X = data$X
#' y = data$class
#'
gaussian_mixture_clusters = function(n = c(50,50), dim = 2, sd.max = .1, sd.noise = .01,
                             r.range = c(.1,1), mixture.range = c(1,4), mixture.sep = .5){
  # k: number of clusters
  # n: number of standpoints

  k = length(n)
  assertthat::assert_that(k>=2)

  # Calculate angles by dividing the circular in k equal parts
  # Calculate centers as points within a sphere of radius in [r.min, r.max]
  Angle = seq(0,2*pi,2*pi/k)[1:k]
  Radius = runif(k, r.range[1], r.range[2])
  Nmixture = sample(mixture.range[1]:mixture.range[2], k, replace = TRUE)

  data = c()
  class = c()
  for (i in 1:k){
    Centers = matrix(0, Nmixture[i], dim)
    Sigma = list()

    theta = Angle[i]
    r = Radius[i]
    for (j in 1:(dim-1))
      Centers[1,j] = r * sin(theta) * cos(theta)^(j-1)
    Centers[1,dim] = r * cos(theta)^(dim-1)

    for (m in 1:Nmixture[i]){
      if (m>1)
        Centers[m,] = Centers[1,] + matrix(runif(dim, -mixture.sep*sd.max, mixture.sep*sd.max),1)

      # Set random covariance matrices
      sigma = matrix(runif(dim*dim, -sd.max/2, sd.max/2), dim)
      sigma[lower.tri(sigma)] = t(sigma)[lower.tri(sigma)]
      diag(sigma) = runif(dim, sd.max/2, sd.max)
      Sigma[[m]] = sigma / m
    }

    n_mixture = divide_interval_int(n[i], Nmixture[i])
    d = gaussian_clusters_with_param(n = n_mixture, center = Centers, sigma =  Sigma)
    data = rbind(data, d[,1:dim])
    class = c(class, rep(i,n[i]))
  }

  data[,1:dim] = data[,1:dim] + matrix(rnorm(sum(n)*dim, 0, sd.noise), sum(n))

  data = list(X = data, class = class)
  return(data)
}

#' Generate multiview clusters from Gaussian distributions with randomly generated parameters
#'
#' @param n vector of number of data points in each cluster
#' The length of \code{n} should be equal to the number of clusters.
#' @param hidden.dim scaler value of dimensions of the hidden state
#' @param observed.dim vector of number of dimensions of the generate clusters.
#' The length of \code{observed.dim} should be equal to the number of clusters.
#' @param sd.max maximum standard deviation of clusters
#' @param sd.noise standard deviation of the added noise
#' @param hidden.r.range the range (min, max) of distance of cluster centers from the origin in the hidden space.
#'
#' @value
#' a list of data points (X) and cluster labels (class)
#'
#' @examples
#' data = multiview_clusters()
#'
multiview_clusters = function (n = c(50,50), hidden.dim = 2, observed.dim = c(2,2,3),
                               sd.max = .1, sd.noise = .01, hidden.r.range = c(.1,1)){

  Nobservation = length(observed.dim)
  k = length(n)
  N = sum(n)

  # Generate hidden states
  hiddenState = gaussian_clusters (n = n, dim = hidden.dim,
                                   sd.max = sd.max, sd.noise = sd.noise, r.range = hidden.r.range)
  class = hiddenState$class
  hiddenState = as.matrix(hiddenState$X)

  # Generate observations with linear trasformation of hidden states
  data_list = list()
  for (i in 1:Nobservation){
    dim_i = observed.dim[i]
    X = matrix(0, N, dim_i)
    W = matrix(runif(hidden.dim*dim_i, .1,1), hidden.dim)
    for (d in 1:dim_i){
      X = hiddenState %*% W
      X = X + matrix(rnorm(N*dim_i, 0, sd.noise), N)
    }
    data_list[[i]] = X
  }

  data = list()
  data[["observation"]] = data_list
  data[["hidden"]] = hiddenState
  data[["class"]] = class
  return(data)
}
