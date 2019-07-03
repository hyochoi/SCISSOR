#' Normalizing log data
#'
#' This function reduces the variation associated with sequentcing depths
#' remaining after a log transformation. The output is the main data object for the downstream
#' shape change detection \link{miscGlobal}.
#'
#' @param inputData log-transformed data matrix (e.g. object from \link{logtransform_data})
#' @param pileup raw coverage data matrix
#' @param exonset exon/intron annotation matrix from \link{process_data}
#' @param smoothness
#' @param makePlot logical, whether to plot variation with respect to overall expression with fitted curve.
#' Default is FALSE
#'
#' @keywords
#' @export
#' @examples
normalize_data = function(inputData,pileup,exonset,
                          smoothness=0.7,
                          makePlot=FALSE, ...) {

  data.centered = center_data(inputData=inputData)
  g1.offset = estimate_offset(data.centered=data.centered,pileup=pileup,
                              exonset=exonset,
                              smoothness=smoothness,makePlot=FALSE)
  msf = data.centered$msf
  datactr = sweep(x=data.centered$outputData,2,g1.offset$g,FUN="/")
  goodcase = g1.offset$goodcase

  ## loop until no different variations
  if (sum((g1.offset$g-1)^2)<1e-10) {
    cat("No adjustment applied. Data do not have enough expression.","\n")
    cat("........Plots are omitted.","\n")
    g2.offset=g1.offset
  } else {
    g2.offset = g1.offset
    k = 1
    while ((k<5) & (max(g2.offset$g)>1.05)) {
      g2.offset = estimate_offset(msf=msf,cenmat=datactr,rawmat=rawmat,exonset=exonset,
                                  smoothness=smoothness,makePlot=F)
      datactr = sweep(x=datactr,2,g2.offset$g,FUN="/")
      k = k+1
    }
    if (makePlot) {
      plot_offset(offset.obj=g1.offset,draw.legend=T,
                  main=paste(GeneName,": before normalization"))
      plot_offset(offset.obj=g2.offset,draw.legend=T,
                  main=paste(GeneName,": after normalization"))
    }

  }
  return(list(outputData=datactr,msf=msf,g1.offset=g1.offset,g2.offset=g2.offset,
              data.center=data.centered$data.center,goodcase=goodcase))
}
