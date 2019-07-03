#' kernel density function
#' @export
kdeplot.hy <- function(x, bandwidth=NULL, indlist=NULL, indlist2=NULL,
                       colmat=NULL,indcol1="red", indcol2="cadetblue",
                       text=FALSE,textwhat=NULL,
                       indpch=16, indcex=0.9, low=0.4, high=0.7, 
                       line.col="black", line.lwd=1,
                       xlab=NULL, ylab=NULL, xlim=NULL, ylim=NULL,main=NULL,cex.main=1.2,
                       ADlegend=FALSE,...){
  ##  Based on 'kdeplotls' from Lingsong's software
  ##  Widen the height of points
  require('KernSmooth') ; require("Lmoments"); require("moments"); require("nortest");
  if (is.null(bandwidth)){
    x.kde <- bkde(x) ;
  }
  if (!is.null(bandwidth)){
    x.kde <- bkde(x, bandwidth=bandwidth) ;
  }
  n <- length(x) ;
  ymin <- quantile(x.kde$y, low)[[1]] ;
  ymax <- quantile(x.kde$y, high)[[1]] ;
  yplot <- ymin + (ymax - ymin)*runif(n) ;
  if (is.null(xlim)){
    xlim=c(min(x), max(x)) ;
  }
  plot(x.kde, type='l', xlim=xlim, ylim=ylim, axes=F, ylab=NA, xlab=NA, col=line.col, lwd=line.lwd) ;
  box() ;
  abline(0, 0, lty=2, col=rgb(.5, .5, .5))
  abline(v=0, lty=2, col=rgb(.5, .5, .5))
  axis(side=1, tck=-.015,labels=NA) ;
  axis(side=1, lwd=0, line=-1, cex=0.2,cex.axis=0.9) ;
  axis(side=2, tck=-0.015,lwd=0,line=-1,cex.axis=0.9)
  # axis(side=2, tck=-.015,labels=NA) ;
  # axis(side=2, lwd=0, line=-1, cex=0.2) ;
  mtext(side=1, xlab, line=1.5, cex=1.3) ;
  mtext(side=2, ylab, line=1.5, cex=1.3) ;
  
  title(main, cex.main=cex.main, font.main=1, line=0.3) ;
  
  if (length(indlist)==0) {indlist=NULL}
  if (is.null(indlist)){
    if (is.null(colmat)){
      colmat <- rep("grey", length(x)) ;
    }
    points(x, yplot, col=colmat, ...) ;
  }
  if (!is.null(indlist)){
    if (is.null(colmat)){
      colmat <- rep("grey", length(x)) ;
      if (!is.null(indlist2)){
        colmat[indlist2] = rep(indcol2,length(indlist2))
      }
      colmat[indlist] <- rep(indcol1, length(indlist)) ;
    }
    indlist0 = unique(indlist,indlist2);
    points(x[-indlist0], yplot[-indlist0], col=colmat[-indlist0], cex=0.8, pch=16) ;
    if (!text){
      if (!is.null(indlist2)){
        points(x[indlist2], yplot[indlist2], col=colmat[indlist2], cex=indcex, pch=indpch) ;
      }
      points(x[indlist], yplot[indlist], col=colmat[indlist], cex=indcex, pch=indpch) ;
    }
    if (text){
      if (!is.null(indlist2)){
        text(x[indlist2], yplot[indlist2], indlist2, col=colmat[indlist2], cex=indcex) ;
      }
      if (is.null(textwhat)) {
        text(x[indlist], yplot[indlist], indlist, col=colmat[indlist], cex=indcex) ;
      } else {
        yplot[indlist] = seq(min(yplot[indlist]),max(yplot[indlist]),length.out=length(indlist))
        text(x[indlist], yplot[indlist], textwhat, col=colmat[indlist], cex=indcex) ;
      }
    }
  }
  if (ADlegend) {
    legend("topleft",bty='n',
           legend=c(paste("skewness =",round(Lmoments(x)[1,3],digits=3)),
                    paste("kurtosis =",round(Lmoments(x)[1,4],digits=3)),
                    paste("A-D      =",round(ADstatWins.hy(x),digits=3))))
  }
}
