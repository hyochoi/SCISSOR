#' Obtain robust mean scale factors
#'
#' This function obtains the mean scale factors
#'
#' @param X an input data matrix
#' @param average
#' @import MASS
#' @export
scale.factor <- function(X,average="mean",trim=0.1,adjval=NULL,robustLM=FALSE) {
  mean.vec = apply(X,1,FUN=function(x){adj.center(x,average=average,trim=trim,adjval=adjval)})
  if (sum(mean.vec^2)>0) {
    if (robustLM) {
      msf = apply(X,2,FUN=function(x){MASS::rlm(x ~ mean.vec - 1,psi=psi.bisquare)$coefficients})
      if (length(which(is.na(msf)==1))>0) {
        temp_id=which(is.na(msf))
        msf[temp_id]=as.vector(t(mean.vec)%*%X[,temp_id]/(sum(mean.vec^2)))
      }
    } else {
      msf = as.vector(t(mean.vec)%*%X/(sum(mean.vec^2)));
    }
  } else {
    msf = rep(0,ncol(X))
  }
  NewX = X - mean.vec%*%t(msf)
  return(list(msf=msf,mean.vec=mean.vec,NewX=NewX))
}
