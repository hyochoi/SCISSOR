#' Plot outlyingness scores with a kernel density estimate
#'
#' @param object For the classes of GSCOutput and LSCOutput, the arguments GSC
#'   and LSC will be ignored.
#' @param GSC logical, whether to plot the OS from GSC. Default is TRUE
#' @param LSC logical, whether to plot the OS from LSC. Default is TRUE
#' @param colmat
#' @param textSC
#'
#' @export
plot_ScissorOS = function(object,GSC=TRUE,LSC=TRUE,
                          colmat=NULL,textSC=TRUE) {
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
    globalCutoff = object$GSCout$cutoff
    localCutoff = object$LSCout$cutoff
  } else {
    if (is(object,"GSCOutput")) {
      n = length(object$OS)
      globalSC = object$SC
      localSC = c()
      globalOS = object$OS
      globalCutoff = object$GSCout$cutoff
      GSC = TRUE
      LSC = FALSE
    } else {
      n = length(object$OS)
      globalSC = object$ignoredCases
      localSC = object$SC
      localOS = object$OS
      localCutoff = object$LSCout$cutoff
      GSC = FALSE
      LSC = TRUE
    }
  }
  if (is.null(colmat)) {
    palette_SCISSOR()
    colmat = rep("darkgrey",n);
    if (length(globalSC)>0) {
      colmat[globalSC] = 1;
    }
    if (length(localSC)>0) {
      colmat[localSC] = 2;
    }
  }
  if (GSC) {
    kdeplot.hy(globalOS,indlist=globalSC,main="Outlyingness Scores from Global",
               text=textSC,high=0.95,low=0.1,colmat=colmat)
    abline(v=globalCutoff)
    legend("topright",bty="n",
           legend=c(paste("# GSCs detected =",length(globalSC))))
  }
  if (LSC) {
    kdeplot.hy(localOS,indlist=localSC,
               main="Outlyingness Scores from Local",
               text=textSC,colmat=colmat,
               high=0.95,low=0.1)
    abline(v=localCutoff)
    legend("topright",bty="n",
           legend=c(paste("# LSCs detected =",length(localSC))))
  }
}
