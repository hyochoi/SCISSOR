#'
#' @export
estimate_offset = function(centerDataResult=NULL,msf=NULL,centeredData=NULL,
                           pileupData,exonset,
                           smoothness=0.7,
                           makePlot=FALSE, ...) {

  if (is.null(centerDataResult)) {
    if ((is.null(msf)) | (is.null(centeredData))) {
      stop("either of msf or centered matrix should be specified.")
    }
  } else {
    msf = centerDataResult$msf;
    centeredData = centerDataResult$outputData;
  }

  n=ncol(centeredData);
  case.sse = apply(centeredData,2,FUN=bisquare.sse)
  exonbase = c()
  for (j in 1:nrow(exonset)) {
    exonbase = c(exonbase,c(exonset[j,2]:exonset[j,3]))
  }
  goodcase = which(apply(pileupData[exonbase,],2,mean)>5)

  ## estimate g by lowess
  if (length(goodcase)>floor(0.1*n)) {
    lowess1 = lowess(case.sse[goodcase] ~ msf[goodcase],f=smoothness)

    tau = lowess1[[2]][which.min((lowess1[[1]]-1)^2)[1]]
    if (tau<0) {
      tau=min(lowess1[[2]][which(lowess1[[2]]>0)])
    }
    ## Rescale case.sse
    case.sse.adj = sqrt(case.sse/tau)
    x = lowess1[[1]][which(lowess1[[2]]>0)]
    y = sqrt(lowess1[[2]][which(lowess1[[2]]>0)]/tau)

    lowess2 = as.list(NULL);
    lowess2$x = msf;
    lowess2$y = interpolate.hy(x=x,y=y,newdata=msf)

    g = lowess2$y;
    g[which(g<1)] = 1;

    if (makePlot) {
      colmat=rep("darkgrey",n)
      colmat[goodcase[order(msf[goodcase],decreasing=TRUE)]] = rainbow(n=length(goodcase),start=0,end=0.8)

      yaxis = yaxis.hy(case.sse.adj)
      plot(msf,case.sse.adj,col=colmat,
           xlab="median scale factor",ylab="variability",
           ylim=c(yaxis[1],yaxis[2]+0.3*(yaxis[2]-yaxis[1])),...)
      lines(lowess2$x[order(lowess2$x)],lowess2$y[order(lowess2$x)],lwd=3,col="black")
      lines(msf[order(msf)],g[order(msf)],lwd=3,col="red")
      abline(v=1,h=1,lty=2)
      # points(msf[order(msf)],g[order(msf)],type="l",lwd=3,col="red")

      legend("topleft",bty="n",
             legend=c(paste("# cases involved =",length(goodcase)),
                      paste("gene length =",length(exonbase)),
                      paste("smoothness  =",round(smoothness,digits=1))))
    }
  } else {
    # No adjustment applied
    g = rep(1,n); case.sse.adj=case.sse;
    lowess2=NULL;

    if (makePlot) {
      colmat=rep("darkgrey",n)
      colmat[goodcase[order(msf[goodcase],decreasing=TRUE)]] = rainbow(n=length(goodcase),start=0,end=0.8)

      yaxis = yaxis.hy(case.sse.adj)
      plot(msf,case.sse.adj,col=colmat,
           xlab="median scale factor",ylab="variability",
           ylim=c(yaxis[1],yaxis[2]+0.3*(yaxis[2]-yaxis[1])), ...)
      abline(v=1,h=1,lty=2)
      # points(msf[order(msf)],g[order(msf)],type="l",lwd=3,col="red")

      legend("topleft",bty="n",
             legend=c(paste("No adjustment applied"),
                      paste("# cases involved =",length(goodcase)),
                      paste("gene length =",length(exonbase))))
    }

  }

  return(list(g=g,msf=msf,sse=case.sse.adj,goodcase=goodcase,
              lowess.fit=lowess2,smoothness=smoothness))
}
# estimate_offset = function(centerDataResult=NULL,msf=NULL,centeredData=NULL,
#                            rawmat,exonset,
#                            smoothness=0.7,makePlot=FALSE,main="Gene") {
#
#   if (is.null(centerDataResult)) {
#     if ((is.null(msf)) | (is.null(centeredData))) {
#       stop("either of msf or centeredData should be specified.")
#     }
#     medscale=msf;
#   } else {
#     medscale = centerDataResult$msf;
#     centeredData = centerDataResult$outdata;
#   }
#
#   n=ncol(centeredData);
#   case.sse = apply(centeredData,2,FUN=bisquare.sse)
#   exonbase = c()
#   for (j in 1:nrow(exonset)) {
#     exonbase = c(exonbase,c(exonset[j,2]:exonset[j,3]))
#   }
#   goodcase = which(apply(rawmat[exonbase,],2,mean)>10)
#
#   ## estimate g by lowess
#   g = rep(1,n)
#   if (length(goodcase)>floor(0.1*n)) {
#     lowess1 = lowess(case.sse[goodcase] ~ medscale[goodcase],f=smoothness)
#
#     tau = lowess1[[2]][which.min((lowess1[[1]]-1)^2)[1]]
#     if (tau<0) {
#       tau=min(lowess1[[2]][which(lowess1[[2]]>0)])
#     }
#     ## Rescale case.sse
#     case.sse.adj = sqrt(case.sse/tau)
#
#     lowess2 = as.list(NULL);
#     lowess2$x = lowess1[[1]][which(lowess1[[2]]>0)]
#     lowess2$y = sqrt(lowess1[[2]][which(lowess1[[2]]>0)]/tau)
#
#     g[goodcase[order(medscale[goodcase])[which(lowess1[[2]]>0)]]] = lowess2$y
#     g[which(g<1)] = 1
#     # g[which(medscale<1)] = 1
#   }
#
#   if (makePlot) {
#     colmat=rep("black",length(goodcase))
#     colmat[order(medscale[goodcase],decreasing=TRUE)] = rainbow(n=length(goodcase),start=0,end=0.8)
#
#     plot(medscale[goodcase],case.sse.adj[goodcase],col=colmat,main=main)
#     lines(lowess2,lwd=3,col="red")
#     abline(v=1,h=1,lty=2)
#     # points(medscale[order(medscale)],g[order(medscale)],type="l",lwd=3,col="red")
#
#     legend("topright",bty="n",
#            legend=c(paste("min =",round(min(lowess2[[2]]),digits=3)),
#                     paste("max =",round(max(lowess2[[2]]),digits=3))))
#     legend("topleft",bty="n",
#            legend=c(paste("# cases involved =",length(goodcase)),
#                     paste("gene length =",length(exonbase)),
#                     paste("smoothness  =",round(smoothness,digits=1))))
#   }
#   return(g)
# }
