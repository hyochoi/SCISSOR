#'
#' @export
compScales = function(x,
                      rmZeroes=FALSE,maxRatio=NULL,precScale=1e-10){
  # Computes the scales sa and sb (above and below the median).
  # Assumes that x is an array of numbers.
  #
  x = x[!is.na(x)] # we always take out NAs
  temp = fastSplitSample(x)
  xa   = temp$xa
  xb   = temp$xb
  med  = temp$med
  sall = scale1StepM((x-med),precScale=precScale)
  if(rmZeroes){ # reduces breakdown value but yields fewer implosions
    xa = xa[xa > precScale]
    xb = xb[xb > precScale]
  }
  sa = scale1StepM(xa,precScale=precScale)
  sb = scale1StepM(xb,precScale=precScale)
  if(!is.null(maxRatio)){
    if(maxRatio < 2) stop("maxRatio must be at least 2")
    sa = min(c(max(sa,sall/maxRatio,na.rm = TRUE),sall*maxRatio),
             na.rm = TRUE)
    sb = min(c(max(sb,sall/maxRatio,na.rm=TRUE),sall*maxRatio),
             na.rm = TRUE)
  }
  return(list(sa=sa,sb=sb,med=med))
}
