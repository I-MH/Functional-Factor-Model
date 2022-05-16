
hatK <- function(evalues){
  # it computes \hat{K} using a ratio-based estimator described on the main paper
  # args:
  #   evalues: a vector of eigenvalues
  # values: a list 
  #   hat_k: index of the best k
  #   ratio: best k
  k_aux <- evalues[-1]/ evalues[-length(evalues)]
  index_k <- which.min(k_aux)
  return(list(hat_k=index_k, ratio=k_aux) )
}
