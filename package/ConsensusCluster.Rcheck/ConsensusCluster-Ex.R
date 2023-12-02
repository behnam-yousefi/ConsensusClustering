pkgname <- "ConsensusCluster"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('ConsensusCluster')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CC_cluster_count")
### * CC_cluster_count

flush(stderr()); flush(stdout())

### Name: CC_cluster_count
### Title: Count the number of clusters based on stability score.
### Aliases: CC_cluster_count

### ** Examples

X = gaussian_clusters()$X
Adj = adj_mat(X, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster=3, max.itter=10)
Result = CC_cluster_count(CM, plot.cdf=FALSE)




cleanEx()
nameEx("Logit")
### * Logit

flush(stderr()); flush(stdout())

### Name: Logit
### Title: Logit function
### Aliases: Logit

### ** Examples

y = Logit(0.5)




cleanEx()
nameEx("adj_conv")
### * adj_conv

flush(stderr()); flush(stdout())

### Name: adj_conv
### Title: Convert adjacency function to the affinity matrix
### Aliases: adj_conv

### ** Examples

Adj_mat = rbind(c(0.0,0.9,0.0),
                c(0.9,0.0,0.2),
                c(0.0,0.2,0.0))
adj_conv(Adj_mat)





cleanEx()
nameEx("adj_mat")
### * adj_mat

flush(stderr()); flush(stdout())

### Name: adj_mat
### Title: Covert data matrix to adjacency matrix
### Aliases: adj_mat

### ** Examples

X = gaussian_clusters()$X
Adj = adj_mat(X, method = "euclidian")




cleanEx()
nameEx("coCluster_matrix")
### * coCluster_matrix

flush(stderr()); flush(stdout())

### Name: coCluster_matrix
### Title: Calculate the Co-cluster matrix for a set of given set of
###   clustering results.
### Aliases: coCluster_matrix

### ** Examples

Clustering = cbind(c(1,1,1,2,2,2),
                   c(1,1,2,1,2,2))
coCluster_matrix(Clustering)




cleanEx()
nameEx("connectivity_matrix")
### * connectivity_matrix

flush(stderr()); flush(stdout())

### Name: connectivity_matrix
### Title: Build connectivity matrix
### Aliases: connectivity_matrix

### ** Examples

connectivity_mat([1,1,1,2,2,2])




cleanEx()
nameEx("consensus_matrix")
### * consensus_matrix

flush(stderr()); flush(stdout())

### Name: consensus_matrix
### Title: Calculate consensus matrix for data perturbation consensus
###   clustering
### Aliases: consensus_matrix

### ** Examples

X = gaussian_clusters()$X
Adj = adj_mat(X, method = "euclidian")
CM = consensus_matrix(Adj, max.cluster=3, max.itter=10, verbos = FALSE)




cleanEx()
nameEx("gaussian_clusters")
### * gaussian_clusters

flush(stderr()); flush(stdout())

### Name: gaussian_clusters
### Title: Generate clusters of data points from Gaussian distribution with
###   randomly generated parameters
### Aliases: gaussian_clusters

### ** Examples

data = gaussian_clusters()
X = data$X
y = data$class




cleanEx()
nameEx("gaussian_clusters_with_param")
### * gaussian_clusters_with_param

flush(stderr()); flush(stdout())

### Name: gaussian_clusters_with_param
### Title: Generate clusters of data points from Gaussian distribution with
###   given parameters
### Aliases: gaussian_clusters_with_param

### ** Examples

center = rbind(c(0,0),
               c(1,1))
sigma = list(diag(c(1,1)),
             diag(2,2))
gaussian_clusters_with_param(c(10, 10), center, sigma)




cleanEx()
nameEx("gaussian_mixture_clusters")
### * gaussian_mixture_clusters

flush(stderr()); flush(stdout())

### Name: gaussian_mixture_clusters
### Title: Generate clusters of data points from Gaussian-mixture-model
###   distributions with randomly generated parameters
### Aliases: gaussian_mixture_clusters

### ** Examples

data = gaussian_mixture_clusters()
X = data$X
y = data$class




cleanEx()
nameEx("generate_gaussian_data")
### * generate_gaussian_data

flush(stderr()); flush(stdout())

### Name: generate_gaussian_data
### Title: Generate a set of data points from Gaussian distribution
### Aliases: generate_gaussian_data

### ** Examples

generate_gaussian_data(10, center=c(0,0), sigma=diag(c(1,1)), label=1)





cleanEx()
nameEx("hir_clust_from_adj_mat")
### * hir_clust_from_adj_mat

flush(stderr()); flush(stdout())

### Name: hir_clust_from_adj_mat
### Title: Hierarchical clustering from adjacency matrix
### Aliases: hir_clust_from_adj_mat

### ** Examples

Adj_mat = rbind(c(0.0,0.9,0.0),
                c(0.9,0.0,0.2),
                c(0.0,0.2,0.0))
hir_clust_from_adj_mat(Adj_mat)




cleanEx()
nameEx("multi_cluster_gen")
### * multi_cluster_gen

flush(stderr()); flush(stdout())

### Name: multi_cluster_gen
### Title: Multiple cluster generation
### Aliases: multi_cluster_gen

### ** Examples

X = gaussian_clusters()$X
cluster_func = function(X, k){return(stats::kmeans(X, k)$cluster)}
Clusters = multi_cluster_gen(X, cluster_func, param = c(2,3))





cleanEx()
nameEx("multi_kmeans_gen")
### * multi_kmeans_gen

flush(stderr()); flush(stdout())

### Name: multi_kmeans_gen
### Title: Multiple K-means generation
### Aliases: multi_kmeans_gen

### ** Examples

X = gaussian_clusters()$X
Clusters = multi_kmeans_gen(X)




cleanEx()
nameEx("multi_pam_gen")
### * multi_pam_gen

flush(stderr()); flush(stdout())

### Name: multi_pam_gen
### Title: Multiple PAM (K-medoids) generation
### Aliases: multi_pam_gen

### ** Examples

X = gaussian_clusters()$X
Clusters = multi_pam_gen(X)




cleanEx()
nameEx("multiview_cluster_gen")
### * multiview_cluster_gen

flush(stderr()); flush(stdout())

### Name: multiview_cluster_gen
### Title: Multiview cluster generation
### Aliases: multiview_cluster_gen

### ** Examples

data = multiview_clusters (n = c(40,40,40), hidden.dim = 2, observed.dim = c(2,2,2),
sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))
X_observation = data[["observation"]]
cluster_func = function(X, k){return(stats::kmeans(X, k)$cluster)}
Clusters = multiview_cluster_gen(X_observation, func = cluster_func)




cleanEx()
nameEx("multiview_clusters")
### * multiview_clusters

flush(stderr()); flush(stdout())

### Name: multiview_clusters
### Title: Generate multiview clusters from Gaussian distributions with
###   randomly generated parameters
### Aliases: multiview_clusters

### ** Examples

data = multiview_clusters()




cleanEx()
nameEx("multiview_consensus_matrix")
### * multiview_consensus_matrix

flush(stderr()); flush(stdout())

### Name: multiview_consensus_matrix
### Title: Calculate consensus matrix for multi-data consensus clustering
### Aliases: multiview_consensus_matrix

### ** Examples

data = multiview_clusters (n = c(40,40,40), hidden.dim = 2, observed.dim = c(2,2,2),
sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))
X_observation = data[["observation"]]
Adj = list()
for (i in 1:length(X_observation))
  Adj[[i]] = adj_mat(X_observation[[i]], method = "euclidian")
CM = multiview_consensus_matrix(Adj, max.cluster = 4, verbos = FALSE)




cleanEx()
nameEx("multiview_kmeans_gen")
### * multiview_kmeans_gen

flush(stderr()); flush(stdout())

### Name: multiview_kmeans_gen
### Title: Multiview K-means generation
### Aliases: multiview_kmeans_gen

### ** Examples

data = multiview_clusters (n = c(40,40,40), hidden.dim = 2, observed.dim = c(2,2,2),
sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))
X_observation = data[["observation"]]
Clusters = multiview_kmeans_gen(X_observation)




cleanEx()
nameEx("multiview_pam_gen")
### * multiview_pam_gen

flush(stderr()); flush(stdout())

### Name: multiview_pam_gen
### Title: Multiview PAM (K-medoids) generation
### Aliases: multiview_pam_gen

### ** Examples

data = multiview_clusters (n = c(40,40,40), hidden.dim = 2, observed.dim = c(2,2,2),
sd.max = .1, sd.noise = 0, hidden.r.range = c(.5,1))
X_observation = data[["observation"]]
Clusters = multiview_pam_gen(X_observation)




cleanEx()
nameEx("pam_clust_from_adj_mat")
### * pam_clust_from_adj_mat

flush(stderr()); flush(stdout())

### Name: pam_clust_from_adj_mat
### Title: PAM (k-medoids) clustering from adjacency matrix
### Aliases: pam_clust_from_adj_mat

### ** Examples

Adj_mat = rbind(c(0.0,0.9,0.0),
                c(0.9,0.0,0.2),
                c(0.0,0.2,0.0))
pam_clust_from_adj_mat(Adj_mat)




cleanEx()
nameEx("spect_clust_from_adj_mat")
### * spect_clust_from_adj_mat

flush(stderr()); flush(stdout())

### Name: spect_clust_from_adj_mat
### Title: Spectral clustering from adjacency matrix
### Aliases: spect_clust_from_adj_mat

### ** Examples

Adj_mat = rbind(c(0.0,0.9,0.0),
                c(0.9,0.0,0.2),
                c(0.0,0.2,0.0))
hir_clust_from_adj_mat(Adj_mat)




### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
