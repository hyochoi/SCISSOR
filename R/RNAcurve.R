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

RNAcurveOne.minfo = function(datmat,indlist,dai,barcode,
                             minfo=NULL,mlegend=FALSE,mcol=NULL,
                             colmat=NULL, same.yaxis=TRUE, is.box=TRUE,
                             title.id=TRUE, title="Sample", title.cex=1,
                             mean.plot=TRUE, mean.col="grey",ylim=NULL,xlim=NULL,
                             cex.axis=0.9,
                             ep.col=NULL,title.barcode=TRUE,
                             ex.num=NULL,int.num=NULL,ex.num.col=NULL,int.num.col=NULL,
                             ylab=NULL, xlab=NULL,yaxis.logcount=NULL,...) {
  ## % Draw individual curves with the mean curve of datmat and vertical lines on exon edges.
  ## % datmat        : d by n  (columns are each sample)
  ## % indlist       : a vector of lists of sample id to be drawn
  ## % mutation.info : mutation information for one gene to be plotted.
  ##                    - all mutations will be indicated by vertical lines
  ##                    - If NULL, no action. (only curve)
  ## % colmat        : color matrix to be used. (length should be same as n)
  ## % ex.num        : exon list that will be colored. (from the right in figure)
  ## % int.num       : intron list that will be colored. (from the right in figure)
  ## % Updated       : September 15, 2017
  ## % Hyo Young Choi

  # Start
  # Extract variables from "gene.input" list.
  epl = dai$epm[,2]; epr = dai$epm[,3]
  num.intron = dai$intron.len;

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
  if (is.null(ex.num.col)){ ex.num.col=candicol3[1] } else {ex.num.col=candicol[ex.num.col]}
  if (is.null(int.num.col)){ int.num.col=candicol3[3] } else {int.num.col=candicol[int.num.col]}


  #  Create colmat by using SVD
  if (is.null(colmat)) {
    x.mda <- svd(datmat) ;
    projmat <- diag(x.mda$d)%*%t(x.mda$v) ;
    projmat[1,] <- -projmat[1,] ;
    colmat <- rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
  }

  # mean.curve = apply(datmat, 1, mean) ;
  mean.curve = apply(datmat, 1, median) ;

  if (title.id) {
    plot.title = paste0(title," # ")
    if (title.barcode){
      plot.title = paste0(plot.title,indlist," | ",barcode[indlist]);
    } else {
      plot.title = paste0(plot.title,indlist);
    }
  } else {
    if (title.barcode) {
      plot.title = paste0(barcode[indlist]);
    } else {
      plot.title = title;
    }
  }

  if (!is.null(yaxis.logcount)) {
    labels = c(0,1,5,7,8,10,20,50,100,500,1000,2000,5000,10000,15000,20000,25000)
    tick.at = log10(labels+yaxis.logcount)-log10(yaxis.logcount);
  } else {
    tick.at = NULL;
    labels = TRUE;
  }

  if (same.yaxis) {
    if (is.null(xlim)) {
      xlim = c(0,nrow(datmat))
    }
    if (is.null(ylim)) {
      ylim0 = yaxis.hy(datmat)
    } else {
      ylim0=ylim
    }
    plot(mean.curve, type='l', lty=2, lwd=0.5, ylim=c(min(0,ylim0[1]),ylim0[2]),
         xlim=xlim, axes=F, ylab=NA, xlab=NA, xaxs="i",yaxs="i", col="white") ;
  } else {
    if (is.null(ylim)){
      ylim=yaxis.hy(datmat[,indlist]);
    }
    if (is.null(xlim)) {
      xlim = c(0,nrow(datmat))
    }
    plot(mean.curve, type='l', lty=2, lwd=0.5, ylim=c(min(0,ylim[1]),ylim[2]),
         xlim=xlim, axes=F, ylab=NA, xlab=NA, xaxs="i",yaxs="i", col="white") ;
  }
  for (i in 1:length(epl)){
    polygon(x=c(rep(epl[i],2),rep(epr[i],2)),y=c(-10000,(max(datmat)+10000),(max(datmat)+10000),-10000),col=ep.col,border=NA) ;
  }
  for (l in ex.num){
    polygon(x=c(rep(epl[l],2),rep(epr[l],2)),y=c(-10000,(max(datmat)+10000),(max(datmat)+10000),-10000),col=ex.num.col,border=NA) ;
  }
  for (l in int.num){
    polygon(x=c(rep(epr[l],2),rep(epl[l+1],2)),y=c(-10000,(max(datmat)+10000),(max(datmat)+10000),-10000),col=int.num.col,border=NA) ;
  }
  abline(v=epl,lty=1,col="lightyellow3",lwd=0.1) ;
  abline(v=epr,lty=1,col="lightyellow3",lwd=0.1) ;
  if (is.box) {
    box()
  }
  title(plot.title, cex.main=title.cex,font.main=1,line=0.3);
  axis(side=1, tck=-0.005, labels=NA) ;
  axis(side=1, lwd=0, line=-0.5, cex.axis=cex.axis) ;
  axis(side=2, tck=-0.015, at=tick.at, labels=labels,lwd=0,line=-1,cex.axis=cex.axis)
  # axis(side=2,lwd=0,line=-0.5,cex.axis=0.9) ;
  # axis(side=2, lwd=0, line=-0.5, cex.axis=0.9) ;
  if (mean.plot){
    points(mean.curve, type='l', lty=1, lwd=2, col=mean.col) ;
  }
  points(datmat[,indlist], type='l', lty=1, col=colmat[indlist], ...) ;
  mtext(side=1, xlab, line=1.5, cex=1.3) ;
  mtext(side=2, ylab, line=1.5, cex=1.3) ;


  if (!is.null(minfo)) {
    mutation.info = minfo[,c(8,9,6)]; mut_start=NULL; mut_end=NULL;
    if (nrow(mutation.info)>0) {
      mut_start = apply(matrix(as.numeric(minfo[,3]),ncol=1),1,
                        FUN=function(x){mutpos.hy(x,exon=exon,is.intron=TRUE,num.intron=dai$intron.len)})
      mut_end = apply(matrix(as.numeric(minfo[,4]),ncol=1),1,
                      FUN=function(x){mutpos.hy(x,exon=exon,is.intron=TRUE,num.intron=dai$intron.len)})
    }
    mutation.info = cbind(mutation.info,mut_start,mut_end);

    case.barcode = barcode[indlist];
    mut_types_temp =  c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
                        "Nonsense_Mutation","Nonstop_Mutation","RNA","Silent","Splice_Site","Translation_Start_Site");
    if (is.null(mcol)) {
      colmut = rep("black",length(mut_types_temp));
      colmut[6] = "blue"; colmut[10] = "red";
    } else {
      colmut = mcol;
    }

    f1=function(x,barcodeset){which(grepl(x,barcodeset)==1)}
    if (nrow(mutation.info)>0) {
      mutype0 = mutation.info[unlist(sapply(case.barcode,f1,mutation.info$barcode2)),3]
      mutpos0 = mutation.info[unlist(sapply(case.barcode,f1,mutation.info$barcode2)),4]
      # mutype0 = mutation.info[which(mutation.info$barcode2==case.barcode),3];
      # mutpos0 = mutation.info[which(mutation.info$barcode2==case.barcode),4];
      if (length(mutpos0)>0) {
        mutpos1 = mutpos0 ;
        colvec = colmut[apply(matrix(mutype0,ncol=1),1,FUN=function(x){which(mut_types_temp==x)})];
        abline(v=mutpos1,lwd=2,col=colvec) ;
        if (mlegend) {
          mutype1 = unique(mutype0); colvec1 = unique(colvec);
          legend("topleft",bty="n",legend=mutype1,lty=1,col=colvec1,cex=1.5)
        }# if (mlegend)
      }
    }# (nrow(mutation.info)>0)
  }
}

##################################################################################
RNAcurve.minfo = function(datmat,indlist=NULL,dai,barcode,
                          mpar=NULL,one.plot=FALSE,same.yaxis=TRUE,is.box=TRUE,
                          minfo=NULL,mlegend=TRUE,mcol=NULL,
                          mzoom=FALSE,mtype.zoom="Splice_Site",
                          colmat=NULL,
                          title.id=TRUE, title="Sample", title.cex=1, title.barcode=TRUE,
                          mean.plot=TRUE, mean.col="grey", fullpar=FALSE,ylim=NULL,xlim=NULL,
                          ep.col=NULL,cex.axis=0.9,
                          ex.num=NULL,int.num=NULL,ex.num.col=NULL,int.num.col=NULL,
                          ylab=NULL, xlab=NULL,yaxis.logcount=NULL,...) {
  ## % Draw individual curves with the mean curve of datmat and vertical lines on exon edges.
  ## % datmat        : d by n  (columns are each sample)
  ## % indlist       : a vector of lists of sample id to be drawn
  ## % mutation.info : mutation information for one gene to be plotted.
  ##                    - all mutations will be indicated by vertical lines
  ##                    - If NULL, no action. (only curve)
  ## % colmat        : color matrix to be used. (length should be same as n)
  ## % ex.num        : exon list that will be colored. (from the right in figure)
  ## % int.num       : intron list that will be colored. (from the right in figure)
  ## % Updated       : September 15, 2017
  ## % Hyo Young Choi

  #  Create colmat by using SVD
  if (is.null(colmat)) {
    x.mda <- svd(datmat) ;
    projmat <- diag(x.mda$d)%*%t(x.mda$v) ;
    projmat[1,] <- -projmat[1,] ;
    colmat <- rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
  }

  epl = dai$epm[,2]; epr = dai$epm[,3]
  num.intron = dai$intron.len;

  if (!is.null(minfo)) {
    mutation.info = minfo[,c(8,9,6)]; mut_start=NULL; mut_end=NULL;
    if (nrow(mutation.info)>0) {
      mut_start = apply(matrix(as.numeric(minfo[,3]),ncol=1),1,
                        FUN=function(x){mutpos.hy(x,exon=exon,is.intron=TRUE,num.intron=dai$intron.len)})
      mut_end = apply(matrix(as.numeric(minfo[,4]),ncol=1),1,
                      FUN=function(x){mutpos.hy(x,exon=exon,is.intron=TRUE,num.intron=dai$intron.len)})
    }
    mutation.info = cbind(mutation.info,mut_start,mut_end);

    # mut_types =  c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation",
    #                "Nonsense_Mutation","Nonstop_Mutation","RNA","Silent","Splice_Site","Translation_Start_Site");
  }

  if (one.plot) {
    if (is.null(indlist)) {
      case=1
      RNAcurveOne.minfo(datmat=datmat,indlist=case,dai,barcode,
                        minfo=NULL,mlegend=FALSE,mcol=mcol,
                        colmat=colmat, same.yaxis=TRUE,is.box=is.box,
                        title.id=FALSE,title=title, title.cex=title.cex,
                        mean.plot=mean.plot, mean.col=,mean.col,ylim=ylim,xlim=xlim,
                        ep.col=ep.col,title.barcode=FALSE,cex.axis=cex.axis,
                        ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                        ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
      for (case in 2:ncol(datmat)) {
        points(datmat[,case],type='l', lty=1, col=colmat[case], ...)
      }
    } else {
      RNAcurveOne.minfo(datmat=datmat,indlist=indlist[1],dai=dai,barcode=barcode,
                        minfo=minfo,mlegend=FALSE,mcol=mcol,
                        colmat=colmat, same.yaxis=TRUE,is.box=is.box,
                        title.id=FALSE,title=title, title.cex=title.cex,
                        mean.plot=mean.plot, mean.col=,mean.col,ylim=ylim,xlim=xlim,
                        ep.col=ep.col,title.barcode=FALSE,cex.axis=cex.axis,
                        ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                        ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
      for (case in indlist[2:length(indlist)]) {
        points(datmat[,case],type='l', lty=1, col=colmat[case], ...)
      }
    }
  } else {
    ncurve = length(indlist);
    if (ncurve==1) {
      if (mzoom) {
        case.barcode = barcode[indlist];
        if (!is.null(minfo)) {
          if (nrow(mutation.info)>0) {
            mutype0 = mutation.info[which(mutation.info$barcode2==case.barcode),3];
            mutpos0 = mutation.info[which(mutation.info$barcode2==case.barcode),4];
            mutpos0 = mutpos0[which(mutype0==mtype.zoom)][1]
            if (length(mutpos0)>0) {
              mutpos1 = mutpos0 ;
              xlim = c(mutpos1-200,mutpos1+200)
            }
          }# (nrow(mutation.info)>0)
        }
      }# if (mzoom)

      RNAcurveOne.minfo(datmat=datmat,indlist=indlist,dai=dai,barcode=barcode,
                        minfo=minfo,mlegend=mlegend,mcol=mcol,
                        colmat=colmat,
                        title.id=title.id,title=title, title.cex=title.cex,
                        mean.plot=mean.plot, mean.col=mean.col,
                        same.yaxis=same.yaxis,ylim=ylim,xlim=xlim,
                        ep.col=ep.col,title.barcode=title.barcode,cex.axis=cex.axis,
                        ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                        ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
    } else {
      if (fullpar){
        par(mfrow=c(6,5), mar=c(2,3,1.5,0.5))
        for (id in indlist){
          if (mzoom) {
            case.barcode = barcode[indlist];
            if (!is.null(minfo)) {
              if (nrow(mutation.info)>0) {
                mutype0 = mutation.info[which(mutation.info$barcode2==case.barcode),3];
                mutpos0 = mutation.info[which(mutation.info$barcode2==case.barcode),4];
                mutpos0 = mutpos0[which(mutype0==mtype.zoom)][1]
                if (length(mutpos0)>0) {
                  mutpos1 = mutpos0 ;
                  xlim = c(mutpos1-200,mutpos1+200)
                }
              }# (nrow(mutation.info)>0)
            }
          }# if (mzoom)

          RNAcurveOne.minfo(datmat=datmat,indlist=id,dai=dai,barcode=barcode,
                            minfo=minfo,mlegend=mlegend,mcol=mcol,
                            colmat=colmat, same.yaxis=same.yaxis,
                            title.id=title.id,title=title, title.cex=title.cex,
                            mean.plot=mean.plot, mean.col=mean.col,ylim=ylim,xlim=xlim,
                            ep.col=ep.col,title.barcode=title.barcode,cex.axis=cex.axis,
                            ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                            ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
        }
      } else {
        if (is.null(mpar)) {
          if ((1 < ncurve) & (ncurve <=4)) {
            par(mfrow=c(1,4), mar=c(2,3,1.5,0.5))
          } else if ((4 < ncurve) & (ncurve <=16)) {
            kk=ceiling(ncurve/4)
            par(mfrow=c(kk,4), mar=c(2,2,1.5,0.5))
          } else if ((16 < ncurve) & (ncurve <=30)) {
            kk=ceiling(ncurve/5)
            par(mfrow=c(kk,5), mar=c(2,3,1.5,0.5))
          } else {
            par(mfrow=c(6,5), mar=c(2,3,1.5,0.5))
          }
        } else {
          par(mfrow=c(6,5), mar=c(2,3,1.5,0.5))
        }
        for (id in indlist){
          if (mzoom) {
            case.barcode = barcode[indlist];
            if (!is.null(minfo)) {
              if (nrow(mutation.info)>0) {
                mutype0 = mutation.info[which(mutation.info$barcode2==case.barcode),3];
                mutpos0 = mutation.info[which(mutation.info$barcode2==case.barcode),4];
                mutpos0 = mutpos0[which(mutype0==mtype.zoom)][1]
                if (length(mutpos0)>0) {
                  mutpos1 = mutpos0 ;
                  xlim = c(mutpos1-200,mutpos1+200)
                }
              }# (nrow(mutation.info)>0)
            }
          }# if (mzoom)

          RNAcurveOne.minfo(datmat=datmat,indlist=id,dai=dai,barcode=barcode,
                            minfo=minfo,mlegend=mlegend,mcol=mcol,
                            colmat=colmat, same.yaxis=same.yaxis,
                            title.id=title.id,title=title, title.cex=title.cex,
                            mean.plot=mean.plot, mean.col=mean.col,ylim=ylim,xlim=xlim,
                            ep.col=ep.col,title.barcode=title.barcode,cex.axis=cex.axis,
                            ex.num=ex.num,int.num=int.num,ex.num.col=ex.num.col,int.num.col=int.num.col,
                            ylab=ylab, xlab=xlab,yaxis.logcount=yaxis.logcount,...)
        }
      }
    }
  }
}

yaxis.hy <- function(mat){
  #  mat : d by n matrix
  tempmax <- max(mat) ;
  tempmin <- min(mat) ;
  templen <- tempmax-tempmin ;
  return(c(tempmin-0.002*templen, tempmax+0.002*templen)) ;
}
