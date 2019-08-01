#' Plot shape changes from SCISSOR
#' @export
plot_ScissorSC = function(object,Data,SCids=NULL,GSC=TRUE,LSC=TRUE,
                          subject.name=NULL,minfo=NULL,colmat=NULL) {

  if (missing(object)) {
    stop("Object is missing. An object for plotting must be specified")
  }
  stopifnot(is(object,c("ScissorOutput","GSCOutput","LSCOutput")))

  if (is(object,"ScissorOutput")) {
    n = ncol(object$logData)
    globalSC = object$globalSC
    localSC = object$localSC
    globalOS = object$GSCout$OS
    localOS = object$LSCout$OS
  } else {
    if (is(object,"GSCOutput")) {
      n = length(object$OS)
      globalSC = object$SC
      localSC = c()
      globalOS = object$OS
      GSC = TRUE
      LSC = FALSE
    } else {
      n = length(object$OS)
      globalSC = object$ignoredCases
      localSC = object$SC
      localOS = object$OS
      GSC = FALSE
      LSC = TRUE
    }
  }

  # If SCids is not NULL, the arguments of GSC and LSC will be ignored.
  if (is.null(SCids)) {
    if ((!GSC) & (!LSC)) {
      stop("No sample was selected for plotting.")
    }
  }
  # Set colors
  if (is.null(colmat)) {
    palette_SCISSOR()
    colmat = rep("black",n);
    if (length(globalSC)>0) {
      colmat[globalSC] = 1;
    }
    if (length(localSC)>0) {
      colmat[localSC] = 2;
    }
  }

  for (case in SCids) {
    RNAcurve.minfo(datmat=object$logData,indlist=case,
                   dai=object$dai,barcode=subject.name,
                   minfo=minfo,
                   lwd=2,colmat=colmat,yaxis.logcount=object$logshiftVal)
  }
}
