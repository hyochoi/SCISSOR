#' written by ------
#'
#' @export
generdir=function(X,ndir){
  # Generates `ndir' directions (with unit norm),
  # orthogonal to d-subsets of X
  # Assumes dim(X) = n x d
  #
  X = data.matrix(X)
  n = dim(X)[1]
  d = dim(X)[2]
  i = 0
  j = 0
  B = array(0.0,dim=c(ndir,d))
  while(i<ndir){
    j = j+1
    #set.seed(j)
    indices = sample(n,d)
    Atemp   = X[indices,,drop=FALSE]
    E       = matrix(1, d, 1)
    if ((qrP = qr(Atemp, tol = 1e-12))$rank == d) {
      B[i+1,] = solve(qrP, E, tol = 1e-12)
      i = i + 1L
    }
  }
  Bnorm  = sqrt(rowSums(B^2))
  Nx     = mean(abs(unname(X)))
  keep   = Bnorm * Nx > 1e-12
  Bnormr = Bnorm[keep]
  B = B[keep, , drop = FALSE]
  A = B/matrix(data = Bnormr, ncol=ncol(B),nrow=nrow(B))
  # rows contain unit norm vectors
  negatives = which(A[,1]<0)
  A[negatives,] = -A[negatives,]
  return(A)
}
