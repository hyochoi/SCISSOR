#'
#' @export
plot_offset = function(offset.obj,draw.legend=TRUE,colmat=NULL,...) {

  g = offset.obj$g; msf = offset.obj$msf; n = length(g);
  goodcase = offset.obj$goodcase;
  sse = offset.obj$sse; lowess.fit = offset.obj$lowess.fit; smoothness=offset.obj$smoothness;

  if (is.null(colmat)) {
    colmat=rep("darkgrey",n)
    colmat[goodcase[order(msf[goodcase],decreasing=TRUE)]] = rainbow(n=length(goodcase),start=0,end=0.8)
  } else {
    colmat[-goodcase] = "darkgrey"
  }

  yaxis = yaxis.hy(sse)
  if (draw.legend) {
    plot(msf,sse,col=colmat,ylim=c(yaxis[1],yaxis[2]+0.2*(yaxis[2]-yaxis[1])),
         xlab="median scale factor",ylab="variation level",...)
  } else {
    plot(msf,sse,col=colmat,ylim=c(yaxis),
         xlab="median scale factor",ylab="variation level",...)
  }
  lines(lowess.fit$x[order(lowess.fit$x)],lowess.fit$y[order(lowess.fit$x)],lwd=2,col="black")
  lines(msf[order(msf)],g[order(msf)],lwd=3,col="red")
  abline(v=1,h=1,lty=2)
  # points(msf[order(msf)],g[order(msf)],type="l",lwd=3,col="red")

  if (draw.legend) {
    legend("topleft",bty="n",
           legend=c(paste("# cases involved =",length(goodcase)),
                    paste("smoothness  =",round(smoothness,digits=1))))
  }

}
