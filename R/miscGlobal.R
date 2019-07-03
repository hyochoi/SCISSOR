#' Detect global shape variants in RNA-seq inputData
#'
#' This function discovers outlying subjects whose RNA-seq have abnormal shapes
#' and provides the most outlying direction for each outlier.
#'
#' @param inputData precessed RNA-seq data from \link{process_data}
#' @param siglev the significance level of the Chi-squared distribution. Default
#'   is 1e-10.
#' @param subt.mean logical, whether to subtract mean before SVD. Default is
#'   FALSE.
#' @param PCnum the number of PCs to be used. If NULL (the default) the number
#'   of PCs will be estimated.
#' @param maxPCnum the maximum number of PCs to be used. Default is 20.
#' @param ADcutoff a cutoff value for checking the normality based on
#'   Anderson-Darling test statistic
#' @param reducedReturn logical, whether to show less results and not to return
#'   large matrices. Default is TRUE.
#'
#' @export
miscGlobal = function(inputData,siglev=1e-5,ADcutoff=3,
                      PCnum=NULL,maxPCnum=20,
                      reducedReturn=TRUE) {

  d = nrow(inputData); n = ncol(inputData);
  z.pca = pca.hy(inputData,subt.mean=FALSE);
  out = pca2misc(scoremat=z.pca$projmat,eigenvalue=z.pca$eigenval,dimension=d,
                 siglev=siglev,ADcutoff=ADcutoff,
                 PCnum=PCnum,maxPCnum=maxPCnum,
                 reducedReturn=FALSE);

  K = out$K;
  if (K==0) {
    residualData = inputData
  } else {
    residualData = inputData - matrix(z.pca$dirmat[,1:K],ncol=K)%*%matrix(z.pca$projmat[1:K,],nrow=K)
  }

  if (K==0) {
    MOD = NULL
    cutoff = 0
  } else if (K==1) {
    MOD = matrix(rep(z.pca$dirmat[,1:K],n),ncol=n);
    cutoff = sqrt(qchisq(p=(1-siglev),df=K))
  } else {
    if (is.null(out$directions)) {
      MOD = NULL;
      cutoff = 0
    } else {
      MOD = z.pca$dirmat[,1:K]%*%out$directions;
      cutoff = sqrt(qchisq(p=(1-siglev),df=K))
    }
  }

  if (reducedReturn) {
    return(list(OS=out$OS,OSpval=out$OSpval,NPS=NULL,MOD=NULL,cutoff=cutoff,
                outliers=out$outliers,outOS=out$outOS,outliers.sort=out$outliers.sort,outOS.sort=out$outOS.sort,
                pca=NULL,residualData=residualData,K=out$K));
  } else {
    return(list(OS=out$OS,OSpval=out$OSpval,NPS=out$NPS,MOD=MOD,cutoff=cutoff,
                outliers=out$outliers,outOS=out$outOS,outliers.sort=out$outliers.sort,outOS.sort=out$outOS.sort,
                pca=z.pca,residualData=residualData,K=out$K));
  }
}
