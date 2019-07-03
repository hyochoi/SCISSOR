#' Plot shape changes from SCISSOR
#' @export
plot_ScissorSC = function(object,SCids=NULL,GSC=TRUE,LSC=TRUE,
                          subject.name=NULL,minfo=NULL,colmat=NULL) {

  # If SCids is not NULL, the arguments of GSC and LSC will be ignored.
  n = ncol(object$datalog)
  outliers1=object$outliers1
  outliers2=object$outliers2
  if (is.null(SCids)) {
    if ((!GSC) & (!LSC)) {
      stop("No sample was selected for plotting.")
    }
    if (GSC) { SCids = c(SCids,outliers1) }
    if (LSC) { SCids = c(SCids,outliers1) }
  }
  # Set colors
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

  for (case in SCids) {
    RNAcurve.minfo(datmat=object$datalog,indlist=case,
                   dai=object$dai,barcode=subject.name,
                   minfo=minfo,
                   lwd=2,colmat=colmat,yaxis.logcount=object$logshift.val)
  }
}
