#' Plot RNA-seq curves
#' @export
RNAcurve = function(datmat,exonset,indlist=NULL,colmat=NULL,
                    plot.title="data",title.cex=1.5,
                    yaxis.logcount=NULL,same.yaxis=TRUE,
                    mean.plot=TRUE,
                    xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,
                    ep.col=NULL,...) {
  if (is.null(colmat)) {
    x.mda <- svd(datmat) ;
    projmat <- diag(x.mda$d)%*%t(x.mda$v) ;
    projmat[1,] <- -projmat[1,] ;
    colmat <- rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
  }
  if (is.null(indlist)) {
    case=1
    RNAcurveOne(datmat=datmat,exonset=exonset,indlist=case,colmat=colmat,
                plot.title=plot.title,title.cex=title.cex,
                mean.plot=mean.plot,yaxis.logcount=yaxis.logcount,
                same.yaxis=same.yaxis,xlim=xlim,ylim=ylim,
                xlab=xlab,ylab=ylab,ep.col=ep.col,...)
    
    for (case in 2:ncol(datmat)) {
      points(datmat[,case],type='l', lty=1, col=colmat[case], ...)
    }
  } else {
    RNAcurveOne(datmat=datmat,exonset=exonset,indlist=indlist[1],colmat=colmat,
                plot.title=plot.title,title.cex=title.cex,
                mean.plot=mean.plot,yaxis.logcount=yaxis.logcount,
                same.yaxis=same.yaxis,xlim=xlim,ylim=ylim,
                xlab=xlab,ylab=ylab,ep.col=ep.col,...)
    for (case in indlist[2:length(indlist)]) {
      points(datmat[,case],type='l', lty=1, col=colmat[case], ...)
    }
  }
  
}
