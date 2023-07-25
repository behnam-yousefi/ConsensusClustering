## Random Data Generation

## 1) Mixture of Gaussians
library(mvtnorm)
library(dplyr)

generateGaussianData <- function(n, center, sigma, label) {
  data = rmvnorm(n, mean = center, sigma = sigma)
  data = data.frame(data)
  names(data) = c("x", "y")
  data = data %>% mutate(class=factor(label))
  data
}

# Two clusters
dataset1 <- {
  # cluster 1
  n = 500
  center = c(5, 5)
  sigma = matrix(c(3, 1, 1, 2), nrow = 2)
  data1 = generateGaussianData(n, center, sigma, 1)
  
  # cluster 2
  n = 500
  center = c(1, 1)
  sigma = matrix(c(1, 0, 0, 1), nrow = 2)
  data2 = generateGaussianData(n, center, sigma, 2)
  
  # all data
  data = bind_rows(data1, data2)
  
  data
}

plot(dataset1[,1:2], pch = 20, col = dataset1[,3], cex = .5)

## N clusters
generateGaussianClusters = function(n, center, sigma){
  # n: vector of number of data points in each cluster (length(n) = Ncluster)
  # center: matrix of centers Ncluster x dim
  # sigma: list of covariance natrices dim X dim, length(sigma) = Ncluster
  
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

n = c(500, 500, 200)
center = rbind(c(5, 5), c(1,1), c(2,4))
sigma = list(matrix(c(3, 1, 1, 2), nrow = 2),
             matrix(c(1, 0, 0, 1), nrow = 2),
             matrix(c(1, 0, 0, 1), nrow = 2))

data = generateGaussianClusters(n = n, center = center, sigma =  sigma)
plot(data[,1:2], pch = 20, col = data[,3], cex = .5)

## 2) Spherical data
library(movMF)

# generate data from Von Mises–Fisher distribution
generateMFData <- function(n, theta, cluster, scale=1) {
  data = rmovMF(n, theta)
  data = data.frame(data[,1], data[,2])
  data = scale * data
  names(data) = c("x", "y")
  data = data %>% mutate(class=factor(cluster))
  data
}

dataset2 <- {
  # cluster 1
  n = 100
  data1 = generateMFData(n, 1 * c(-1, 1), 1, 5)
  
  # cluster 2
  n = 100
  data2 = generateMFData(n, 1 * c(1, -1), 2, 1)
  
  # all data
  data = bind_rows(data1, data2)
  
  data
}

plot(dataset2[,1:2], pch = 20, col = dataset2[,3], cex = .5)

## 3) Spiral data

generateSpiralData <- function(n) {
  maxRadius = 7
  xShift = 2.5
  yShift = 2.5
  angleStart = 2.5 * pi
  noiseVariance = 0.1
  
  # first spiral
  firstSpiral <- function() {
    d1 = data.frame(0:(n-1))
    colnames(d1) <- c("i")
    d1 %>% mutate(angle = angleStart + 2.5 * pi * (i / n),
                  radius = maxRadius * (n + n/5 - i)/ (n + n/5),
                  x = radius * sin(angle),
                  y = radius * cos(angle),
                  class=1)
  }
  d1 = firstSpiral()
  
  # second spiral
  d2 = d1 %>% mutate(x = -x, y = -y, class=2)
  
  # combine, add noise, and shift
  generateNoise <- function(n) {
    sigma = matrix(c(noiseVariance, 0, 0, noiseVariance), nrow = 2)
    noise = rmvnorm(n, mean = c(0, 0), sigma = sigma)
    df = data.frame(noise)
    colnames(df) <- c("xNoise", "yNoise")
    df
  }
  d1 %>%
    bind_rows(d2) %>%
    bind_cols(generateNoise(2*n)) %>%
    transmute(x = x + xShift + xNoise,
              y = y + yShift + yNoise,
              class = factor(class))
}

# dataset3 = generateSpiralData(100)
# plot(dataset3[,1:2], pch = 20, col = dataset3[,3], cex = .5)

