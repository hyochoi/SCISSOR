#' Transforming data
#'
#' This function tansforms the input data
#' @param inputData data to be transformed
#' @param logshiftVal
#' @keywords log shift parameter
#' @export
#' @examples

logtransform_data = function(inputData, logshiftVal=NULL, param.grid=NULL, draw.plot=FALSE) {
  # Find exon bases out of intron-contained coverage
  #     & Take log-transformation.
  if (is.null(logshiftVal)) {
    logshiftVal = getShiftParam(X=inputData, param.grid=param.grid, draw.plot=draw.plot)$optim.param
  }
  outputData = log10(data + logshiftVal) - log10(logshiftVal) ;
  return(list(outputData=outputData,logshiftVal=logshiftVal))
}
