#'
#' @export
interpolate.hy = function(x,y,newdata) {

  inner.fn = function(x,y,newx) {
    x.sorted = x[order(x)]
    y.sorted = y[order(x)]

    # find two values that are closest to newdata and cover newdata between them
    x.diff = x.sorted-newx
    if (max(x.diff)<0) {
      # this case is when newx is beyond the maximum x
      newy = y.sorted[length(y.sorted)]
    } else if (min(x.diff)>0) {
      # this case is when newx is beyond the minimum x
      newy = y.sorted[1]
    } else {
      i = which.min(x.diff^2)
      if (x.diff[i]==0) {
        newy = y.sorted[i]
      } else {
        if (x.diff[i]<0) {
          il = i; iu = i+1;
        } else {
          il = i-1; iu = i;
        }
        if (abs(x.sorted[iu]-x.sorted[il])>1e-5) {
          b = (y.sorted[iu]-y.sorted[il])/(x.sorted[iu]-x.sorted[il])
          newy = y.sorted[il]+b*(newx-x.sorted[il])
        } else {
          newy = y.sorted[il]
        }
      }
    }
    return(newy)
  }
  return(apply(matrix(newdata,ncol=1),1,FUN=function(z){inner.fn(x=x,y=y,newx=z)}))
}

