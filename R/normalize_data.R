#' Normalizing log data
#'
#' This function reduces the variation associated with sequencing depths
#' remaining after a log transformation. The output is the main data object for the downstream
#' shape change detection \link{miscGlobal}.
#'
#' @param inputData log-transformed data matrix (e.g. object from \link{logtransform_data})
#' @param pileupData raw coverage data matrix
#' @param exonset exon/intron annotation matrix from \link{process_data}
#' @param smoothness
#' @param makePlot logical, whether to plot variation with respect to overall expression with fitted curve.
#' Default is FALSE
#' @param Gene gene name to print in the plot
#'
#' @export
normalize_data = function(inputData,
                          pileupData,Ranges,
                          smoothness=0.7,
                          makePlot=FALSE, ...) {

  centerDataResult = center_data(inputData=inputData,Ranges=Ranges)
  exonset = Ranges$lRanges
  Gene = Ranges$Gene
  g1.offset = estimate_offset(centerDataResult=centerDataResult,
                              msf=NULL,centeredData=NULL,
                              pileupData=pileupData,
                              exonset=exonset,
                              smoothness=smoothness,makePlot=FALSE)
  msf = centerDataResult$msf
  normalizedData = sweep(x=centerDataResult$outputData,2,g1.offset$g,FUN="/")
  goodcase = g1.offset$goodcase
  if (makePlot) {
    if (is.null(Gene)) {
      plot_offset(offset.obj=g1.offset,draw.legend=T,
                  main="Before normalization",...)
    } else {
      plot_offset(offset.obj=g1.offset,draw.legend=T,
                  main=paste0(Gene," | Before normalization"),...)
    }
  }

  ## loop until no different variations
  if (sum((g1.offset$g-1)^2)<1e-10) {
    cat("No adjustment applied. Data do not have enough expression.","\n")
    cat("........Plots are omitted.","\n")
    g2.offset=g1.offset
  } else {
    g2.offset = g1.offset
    k = 1
    while ((k<5) & (max(g2.offset$g)>1.05)) {
      g2.offset = estimate_offset(centerDataResult=NULL,
                                  msf=msf,centeredData=normalizedData,
                                  pileupData=pileupData,
                                  exonset=exonset,
                                  smoothness=smoothness,makePlot=FALSE)
      normalizedData = sweep(x=normalizedData,2,g2.offset$g,FUN="/")
      k = k+1
    }
    if (makePlot) {
      if (is.null(Gene)) {
        plot_offset(offset.obj=g2.offset,draw.legend=T,
                    main="After normalization",...)
      } else {
        plot_offset(offset.obj=g2.offset,draw.legend=T,
                    main=paste0(Gene," | After normalization"),...)
      }
    }

  }
  return(list(outputData=normalizedData,msf=msf,
              g1.offset=g1.offset,g2.offset=g2.offset,
              goodcase=goodcase))
}
