#'
#' @export
pca.hy <- function(data, subt.mean=TRUE){
  ##  PCA with the covariance matrix (X%*%t(X)/n-1) where X=data
  ##  data: d by r (d is the number of variables, n is the number of subjects)
  ##  If subt.mean=TRUE, subtract mean from the data.
  if (subt.mean){
    x <- data-apply(data,1,mean) ;
    div.fac <- ncol(x)-1;
  } else {
    x <- data ;
    div.fac <- ncol(x);
  }

  x.mda <- svd(x) ;
  dirmat <- x.mda$u ;
  if (dirmat[1,1]<0){  dirmat[,1] <- -dirmat[,1]}
  projmat <- t(dirmat)%*%x ;
  eigenval <- ((x.mda$d)^2)/div.fac ;
  stdprojmat <- diag(1/sqrt(eigenval))%*%projmat		#  standardized
  return(list(dirmat=dirmat, projmat=projmat, eigenval=eigenval, stdprojmat=stdprojmat))
}
