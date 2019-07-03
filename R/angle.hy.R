#'
#' @export
angle.hy = function(x,y,perp=TRUE){
  
  angle2.hy=function(xx,y,perp=TRUE){
    sumx2 = sum(xx^2); sumy2 = sum(y^2)
    if ((sumx2>0) & (sumy2>0)) {
      ang0 = min(1,as.vector((t(xx)%*%y)/sqrt(sumx2*sumy2)));
    } else {
      ang0 = 0
    }
    ang = acos(ang0)*(180/pi);
    if (perp){
      ang = min(ang,180-ang)
    }
    return(ang)
  }
  X = as.matrix(x); y = as.vector(y);
  angle = apply(X,2,FUN=function(x){angle2.hy(x,y,perp=perp)});
  return(angle);
}
