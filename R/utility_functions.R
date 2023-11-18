# Utilities

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

