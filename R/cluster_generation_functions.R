## Random data generation functions

# library(mvtnorm)
source("R/mvnorm.R")
library(dplyr)


divide_interval_int = function(x,n){
  y = rep(floor(x/n), n-1)
  y[n] = x - sum(y)
  return(y)
}

generateGaussianData = function(n, center, sigma, label) {
  data = rmvnorm(n, mean = center, sigma = sigma)
  data = data.frame(data)
  data = data %>% mutate(class=factor(label))
  return(data)
}

gaussian_clusters_with_param = function(n, center, sigma){
  # n: vector of number of data points in each cluster (length(n) = Ncluster)
  # center: matrix of centers Ncluster x dim
  # sigma: list of covariance matrices dim X dim, length(sigma) = Ncluster
  
  Ncluster = length(n)
  dim = ncol(center)
  assertthat::assert_that(Ncluster == nrow(center))
  assertthat::assert_that(Ncluster == length(sigma))
  
  data = c()
  for (i in 1:Ncluster){
    n_i = n[i]
    center_i = center[i,]
    sigma_i = sigma[[i]]
    data_i = generateGaussianData(n_i, center_i, sigma_i, i)
    data = rbind(data, data_i)
  }
  return(data)
}

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