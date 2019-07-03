#'
#' @export
pos.sd = function(x,x0=NULL,cutoff=0) {
  # Find non-zero entries of x0 
  #   and compute variance of x on those entries
  if (is.null(x0)) {
    y = x[x>cutoff]
  } else {
    y = x[x0>cutoff]
  }
  if (length(y)>0) {
    return(sd(y))
  } else {
    return(0)
  }
}
