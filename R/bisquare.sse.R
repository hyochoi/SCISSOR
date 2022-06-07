#'
bisquare.sse = function(x) {
  sig0 = sd(x)
  if (sig0<1e-10) {
    return(0)
  } else {
    k = 4.685*sig0
    c = (k^2)/6
    bx = rep(c, length(x))
    bx[which(abs(x) <= k)] = c*(1-((1-((bx[which(abs(x) <= k)]/k)^2))^3))
    return(sum(bx))
  }
}
# bisquare.sse = function(x) {
#   sig0 = sd(x)
#   if (sig0<1e-10){
#     return(0)
#   } else {
#     k = 4.685*sig0
#     bisquare.fnt = function(y,k) {
#       c = (k^2)/6
#       if (abs(y)<=k) {
#         return(c*(1-((1-((y/k)^2))^3)))
#       } else {
#         return(c)
#       }
#     }
#     bx = sapply(x,FUN=function(x){bisquare.fnt(y=x,k=k)})
#     return(sum(bx))
#   }
# }

