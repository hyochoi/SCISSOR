#' Determining the most outlying structure using sparse basis for one subject in a data set
#'
#' This function discovers outlying structure for each subject as well as outlyingness scores
#'
#' @param basis a matrix with sparse basis in columns
#' @param data normalized data through pre-processing procedure
#' @param case a subject ID number in the basis matrix
#'
#' @export
get_sparsePD4one = function(basis,data,case) {
  # basis = a matrix with bases in columns (d by k)
  # data = a data matrix with samples in columns (d by n)
  ## % inner functions
  normalize_vector=function(x) {
    return(x/sqrt(sum(x^2)))
  }
  if (is.null(dim(basis))) {
    basis=matrix(basis,ncol=1)
    NPS=pd.rate.hy(as.vector(t(basis)%*%data))
    basis.coef=1
    directions=basis
    if (NPS[case]<0) {
      NPS=-NPS
      basis.coef=-basis.coef
      directions=-directions
    }
  } else {
    n=ncol(data); d=nrow(data); k=ncol(basis);
    x=t(data)%*%basis  # n by k
    ndir=300*k
    A0=generdir(x,ndir=ndir)   # ndir by k
    A1=t(apply((t(basis)%*%data),2,normalize_vector))
    A2=diag(rep(1,k))
    A=rbind(A0,A1,A2)

    ## % calculate OS
    Y = x %*% t(A) # project x onto A (n by ndir)
    ADstat = apply(Y,2,ADstatWins.hy);
    indir = which(ADstat<ADcutoff);
    A = A[indir,];
    Y = x %*% t(A) # project x onto A (n by length(indir))

    out_temp = apply(X=Y,MARGIN=2,pd.rate.hy); # n by length(indir)
    indexmax=min(which.max(abs(out_temp[case,])))
    NPS = out_temp[,indexmax];  # n by n
    basis.coef = t(A[indexmax,])
    directions=basis%*%matrix(basis.coef,ncol=1);
    if (NPS[case]<0) {
      NPS=-NPS
      basis.coef=-basis.coef
      directions=-directions
    }
  }
  return(list(NPS=NPS,directions=directions,basis.coef=basis.coef))
}
