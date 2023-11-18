#' Covert data matrix to adjacency matrix
#'
#' @param x a list obtained by the \code{mnda_embedding_2layer()} function.
#' @param method method for adjusting p-value (including methods on \code{p.adjust.methods}).
#' If set to "none" (default), no adjustment will be performed.
#'
#' @value
#' calculated adjacency matrix from the data matrix using the specified methods
#'
#' @examples
#' X = gaussian_clusters()$X
#' Adj = adj_mat(X, method = "euclidian")
#'
adj_mat = function(x, method = "euclidian"){

  if (method == "cosine"){
    X = X / apply(X, 1, function(x){sqrt(sum(x^2))})
    AdjMat = X %*% t(X)
  }else{
    dists = as.matrix(dist(X, method))
    AdjMat = (max(dists) - dists)/max(dists)
  }
  return(AdjMat)
}

