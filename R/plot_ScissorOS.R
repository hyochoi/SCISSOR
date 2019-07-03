#' Plot outlyingness scores with a kernel density estimate
#' @export
plot_ScissorOS = function(object,GSC=TRUE,LSC=TRUE,
                          colmat=NULL,textSC=TRUE) {

  n = ncol(object$datalog)
  outliers1=object$outliers1
  outliers2=object$outliers2
  if (is.null(colmat)) {
    palette.hy()
    colmat = rep("darkgrey",n);
    if (length(outliers1)>0) {
      colmat[outliers1] = 1;
    }
    if (length(outliers2)>0) {
      colmat[outliers2] = 2;
    }
  }
  if (GSC) {
    kdeplot.hy(object$GSCout$OS,indlist=outliers1,main="Outlyingness Scores from Global",
               text=textSC,high=0.95,low=0.1,colmat=colmat)
    abline(v=object$GSCout$cutoff)
    legend("topright",bty="n",
           legend=c(paste("# GSCs detected =",length(outliers1))))
  }
  if (LSC) {
    kdeplot.hy(object$LSCout$OS,indlist=outliers2,
               main="Outlyingness Scores from Local",
               text=textSC,colmat=colmat,
               high=0.95,low=0.1)
    abline(v=object$LSCout$cutoff)
    legend("topright",bty="n",
           legend=c(paste("# LSCs detected =",length(outliers2)),
                    paste("K-S df     =",object$LSCout$ks.df),
                    paste("K-S pval   =",round(object$LSCout$ks.pval,5))))
  }
}
