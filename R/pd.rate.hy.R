#' Projection depth
#' @export
pd.rate.hy = function(x,qrsc=FALSE) {
  # projection depth
  m = mad(x)
  if (m<1e-5) {
    return(rep(0,length(x)))
  } else {
    if (qrsc) {
      rsc=compScales(x)
      y=rep(0,length(x))

      above.ind=which((x-rsc$med)>=0)
      below.ind=which((x-rsc$med)<0)
      if (rsc$sa>1e-5) {
        y[above.ind]=(x[above.ind]-rsc$med)/rsc$sa
      }
      if (rsc$sb>1e-5) {
        y[below.ind]=(x[below.ind]-rsc$med)/rsc$sb
      }
      return(y)
    } else {
      return((x-median(x))/m)
    }
  }
}
