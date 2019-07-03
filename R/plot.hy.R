#'
#' @export
plot.hy = function(x,y,indlist=NULL,text=F,
                   colmat=NULL,indcol="red",
                   cex=1,indcex=1.2,...) {
  n=length(x)
  if (length(indlist)==0) {indlist=NULL}
  if (is.null(indlist)) {
    if (is.null(colmat)) {
      colmat=rep("black",length(x))
    }
    plot(x=x,y=y,col=colmat,...)
  } else {
    if (is.null(colmat)) {
      colmat=rep("grey",n)
      colmat[indlist]=rep(indcol,length(indlist))
    }
    plot(x[-indlist],y=y[-indlist],col=colmat[-indlist],cex=cex,xlim=yaxis.hy(x),ylim=yaxis.hy(y),...)
    if (!text) {
      points(x[indlist],y[indlist],col=colmat[indlist],cex=indcex,...)
    } else {
      text(x[indlist],y[indlist],indlist,col=colmat[indlist],cex=indcex,...)
    }
  }
}
