RNAcurveOne = function(datmat,exonset,indlist,colmat=NULL,
                       plot.title="data",title.cex=1.5,
                       yaxis.logcount=NULL,same.yaxis=TRUE,
                       mean.plot=TRUE,
                       xlab=NULL,ylab=NULL,xlim=NULL,ylim=NULL,
                       ep.col=NULL,...) {
  
  require("RColorBrewer");
  # display.brewer.all()
  candicol1=c(brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
              brewer.pal(8,"Pastel2")[6],
              brewer.pal(9,"YlOrBr")[1],
              brewer.pal(9,"YlOrRd")[1],
              brewer.pal(9,"YlOrRd")[2],
              brewer.pal(9,"Reds")[1],
              brewer.pal(9,"RdPu")[1],
              brewer.pal(9,"OrRd")[1],
              brewer.pal(9,"OrRd")[2],
              brewer.pal(9,"Oranges")[1],
              brewer.pal(9,"Oranges")[2],
              "aliceblue");
  candicol2=brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
  candicol3=brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
  candicol = c(candicol2,candicol3);
  
  if (is.null(ep.col)){ ep.col = candicol1[9] } else {ep.col=candicol1[ep.col]}
  
  epl = exonset[,2]; epr = exonset[,3];
  
  if (is.null(colmat)) {
    x.mda <- svd(datmat) ;
    projmat <- diag(x.mda$d)%*%t(x.mda$v) ;
    projmat[1,] <- -projmat[1,] ;
    colmat <- rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
  }
  
  if (!is.null(yaxis.logcount)) {
    labels = c(0,1,5,7,8,10,20,50,100,500,1000,2000,5000,10000,15000,20000,25000)
    tick.at = log10(labels+yaxis.logcount)-log10(yaxis.logcount);
  } else {
    tick.at = NULL;
    labels = TRUE;
  }
  
  median.curve = apply(datmat, 1, median) ;
  if (same.yaxis) {
    ylim0 = yaxis.hy(datmat)
    if (is.null(xlim)) {
      xlim = c(0,nrow(datmat))
    } 
    plot(median.curve, type='l', lty=2, lwd=0.5, ylim=c(min(0,ylim0[1]),ylim0[2]), xlim=xlim, xaxs="i",
         axes=F, ylab=NA, xlab=NA, col="white") ;
  } else {
    if (is.null(ylim)){
      ylim=yaxis.hy(datmat[,indlist]);
    }
    if (is.null(xlim)) {
      xlim = c(0,nrow(datmat))
    } 
    plot(median.curve, type='l', lty=2, lwd=0.5, ylim=ylim,  xlim=xlim, xaxs="i",
         axes=F, ylab=NA, xlab=NA, col="white") ;
  }
  
  for (i in 1:length(epl)){
    polygon(x=c(rep(epl[i],2),rep(epr[i],2)),y=c(-10000,(max(datmat)+10000),(max(datmat)+10000),-10000),col=ep.col,border=NA) ;
  }
  box() ;
  
  title(plot.title, cex.main=title.cex,font.main=1,line=0.3);
  abline(h=0, lty=2) ;
  axis(side=1, tck=-0.005, labels=NA) ;
  axis(side=1, lwd=0, line=-1, cex.axis=0.9) ;
  axis(side=2, tck=-0.015, at=tick.at, labels=labels,lwd=0,line=-1,cex.axis=0.9)
  # axis(side=2,lwd=0,line=-0.5,cex.axis=0.9) ;
  # axis(side=2, lwd=0, line=-0.5, cex.axis=0.9) ;
  if (mean.plot){
    points(median.curve, type='l', lty=1, lwd=2, col="grey") ;
  }
  points(datmat[,indlist], type='l', lty=1, col=colmat[indlist], ...) ;
  mtext(side=1, xlab, line=1.5, cex=1.3) ;
  mtext(side=2, ylab, line=1.5, cex=1.3) ;
  
}
