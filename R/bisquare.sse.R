#'
#' @export
bisquare.sse = function(x) {
  sig0 = sd(x)
  if (sig0<1e-10){
    return(0)
  } else {
    k = 4.685*sig0
    bisquare.fnt = function(y,k) {
      c = (k^2)/6
      if (abs(y)<=k) {
        return(c*(1-((1-((y/k)^2))^3)))
      } else {
        return(c)
      }
    }
    bx = apply(matrix(x,ncol=1),1,FUN=function(x){bisquare.fnt(y=x,k=k)})
    return(sum(bx))
  }
}
