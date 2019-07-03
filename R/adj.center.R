#' Provide a center value adjusted by a lower bound
#'
#' This function calculates a center of the values above a given bound.
#'
#' @param x an input vector, data values.
#' @param average type of the center being calculated, with choices "median" and "mean".
#' Default is "median".
#' @param adjval the lower bound above which the values will be included to calculate the center.
#'
#' @export
adj.center = function(x,average="mean",trim=0.1,adjval=NULL) {
  ##  % Obtain adjusted center (mean or median)
  if (is.null(adjval)) {
    if (average=="median") {
      center.val = median(x);
    } else {
      center.val = mean(x,trim=trim);
    }
  } else {
    nonzero.pos = which(x>adjval)
    if (length(nonzero.pos)>0) {
      if (average=="median") {
        center.val = median(x[nonzero.pos])
      } else {
        center.val = mean(x[nonzero.pos])
      }
    } else {
      center.val = 0
    }
  }
  return(center.val);
}
