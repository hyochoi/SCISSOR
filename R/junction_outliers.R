#'
#' @export
find_junction_outliers_i = function(i,JSR.matrix,JSR.table,splice.ratio,
                                    cutoff.a=0.1,cutoff.b=0.9) {

  junctions.g = matrix(as.numeric(split_junction(JSR.table[,1])),ncol=2,byrow=T)
  junction.start = unique(junctions.g[,1])

  index.temp=which(junctions.g[,1]==junctions.g[i,1])
  junction.sum = apply(matrix(JSR.matrix[index.temp,],ncol=n),2,sum);
  # Find outlier candidates
  a0=boxplot.hy(splice.ratio[i,],robust=T,print.plot=F,
                ylim=c(0,max(splice.ratio)),
                title=junctions.g[i],title.cex=0.9,
                indlist=NULL,indlist.col="red",
                outliers.range=iqr.range,outliers.col="grey")
  # Outliers with high evidence
  outliers.tmp1 = outliers.tmp2 = c()
  if ((length(which(splice.ratio[i,]<cutoff.b))<0.1*n) & (length(a0$outliers.below)>0)) {
    ## This is the case of possible exon skipping / cryptic events / SV
    outliers.tmp1 = a0$outliers.below[which((splice.ratio[i,a0$outliers.below]-ES_cutoff_fn(junction.sum[a0$outliers.below]))<0)]
    # outliers.tmp1 = outliers.tmp1[sapply(outliers.tmp1,FUN=function(k) {(JSR.matrix[i,k]<(0.8*median(JSR.matrix[which(JSR.matrix[,k]>5),k])))})]
  }
  if ((length(which(splice.ratio[i,]>cutoff.a))<0.1*n) & (length(a0$outliers.above)>0)) {
    ## This is the case of possible alternative exon / cryptic events / SV
    outliers.tmp2 = a0$outliers.above[which(splice.ratio[i,a0$outliers.above]>=cutoff.a)]
    outliers.tmp2 = outliers.tmp2[which(sapply(outliers.tmp2,FUN=function(k) { (JSR.matrix[i,k]>max(median(JSR.matrix[which(JSR.matrix[,k]>=5),k])*0.2,10))})==T)]
  }

  output = rbind(matrix(cbind(outliers.tmp1,rep("-",times=length(outliers.tmp1)),matrix(JSR.table[rep(i, each = length(outliers.tmp1)), ],nrow=length(outliers.tmp1))),ncol=(2+ncol(JSR.table))),
                 matrix(cbind(outliers.tmp2,rep("+",times=length(outliers.tmp2)),matrix(JSR.table[rep(i, each = length(outliers.tmp2)), ],nrow=length(outliers.tmp2))),ncol=(2+ncol(JSR.table))))
  if (dim(output)[1]>0) {
    colnames(output) = c("outliers","sign",colnames(JSR.table))
    rownames(output) = rep(rownames(JSR.table)[i],times=dim(output)[1])
  } else {
    output = NULL
  }
  return(output)
}

#'
#' @export
ES_cutoff_fn = function(x) {
  ## x = c(1,5,10,20,50,100,500,1000)
  ## cbind(x,ES_cutoff_fn(x))
  tmp_fn = function(t,a,b) {
    slope = (b[2]-a[2])/(log10(b[1])-log10(a[1]))
    ycoef = slope*log10(a[1]) - a[2]
    return(slope*t-ycoef)
  }
  inner_fn = function(y) {
    if (y<10) {
      return(0)
    } else if (y<=20) {
      return(tmp_fn(t=log10(y),a=c(10,0.2),b=c(20,0.6)))
    } else if (y<=50) {
      return(tmp_fn(t=log10(y),a=c(20,0.6),b=c(50,0.7)))
    } else if (y<=100) {
      return(tmp_fn(t=log10(y),a=c(50,0.7),b=c(100,0.8)))
    } else if (y<=500) {
      return(tmp_fn(t=log10(y),a=c(100,0.8),b=c(500,0.9)))
    } else if (y<=1000) {
      return(tmp_fn(t=log10(y),a=c(500,0.9),b=c(1000,0.95)))
    } else {
      return(0.95)
    }
  }
  return(sapply(x,inner_fn))
}

#'
#' @export
boxplot.hy=function(value,robust=FALSE,
                    title=NULL,title.cex=2,
                    ylab=NULL,
                    value.labels=NULL,tick.at=NULL,
                    ylim=NULL,
                    indlist=NULL,indlist.col=NULL,
                    indcex=2,text=FALSE,textwhat=NULL,
                    print.outliers=FALSE,outliers.range=NULL,outliers.text=FALSE,
                    colmat=NULL,outliers.col=NULL,
                    box.col="black",print.plot=TRUE) {
  require(wesanderson); require(RColorBrewer)
  n = length(value)
  qvalue = quantile(value)
  iqr = qvalue[4]-qvalue[2]
  if (!robust) {
    if (is.null(outliers.range)) {
      fence=c(qvalue[2]-iqr*1.5,qvalue[4]+iqr*1.5)
    } else {
      fence=c(qvalue[2]-iqr*outliers.range,qvalue[4]+iqr*outliers.range)
    }
  } else {
    iqr.u=qvalue[4]-qvalue[3] # iqr for large values
    iqr.b=qvalue[3]-qvalue[2] # iqr for low values
    if (is.null(outliers.range)) {
      fence=c(qvalue[2]-2*iqr.b*1.5,qvalue[4]+2*iqr.u*1.5)
    } else {
      fence=c(qvalue[2]-2*iqr.b*outliers.range,qvalue[4]+2*iqr.u*outliers.range)
    }
  }
  box.outliers.below=which(value<fence[1])
  box.outliers.above=which(value>fence[2])
  box.outliers = c(which(value<fence[1]),which(value>fence[2]))

  if (print.outliers) {
    # indlist=box.outliers
    if (outliers.text) {
      text=TRUE
    }
  }
  if (print.plot) {
    if ((!is.null(indlist)) & (length(indlist)==0)) {indlist=NULL}
    # Define colors for points
    if (is.null(colmat)) {
      colmat=rep("grey",n)
      if (!is.null(box.outliers)) {
        if (is.null(outliers.col)) {
          outliers.col = "red"
        }
        colmat[box.outliers] = outliers.col;
      }

      if (!is.null(indlist)) {
        if (is.null(indlist.col)) {
          indlist.col=wes_palette(n=5,"Darjeeling1")[2]
        }
        colmat[indlist] = indlist.col;
      }
    }

    # Define x values
    x = rnorm(n=length(value),mean=0.5,sd=0.03)

    # Define ylim if it is null
    if (is.null(ylim)) {
      ylim=yaxis.hy(value)
      if (ylim[1]==ylim[2]) {
        if (ylim[1]==0) {
          ylim[1]=-0.5; ylim[1]=0.5
        } else {
          ylim[1]=ylim[1]*0.5; ylim[2]=ylim[2]*1.5
        }
      }
    }
    if (is.null(indlist)) {
      plot(x,value,axes=F,col=colmat,ylim=ylim,
           xlab=NA,ylab=NA,xlim=c(min(x)-0.05,max(x)+0.05),pch=19,cex=1,lwd=3)
    } else {
      plot(x[-indlist],value[-indlist],axes=F,col=colmat[-indlist],ylim=ylim,
           xlab=NA,ylab=NA,xlim=c(min(x)-0.05,max(x)+0.05),pch=19,cex=1,lwd=3)
    }

    if (is.null(value.labels)) {
      tick.at=seq(ylim[1],ylim[2],length=5)
      if ((qvalue[5]-qvalue[1]>0) & (qvalue[5]<1e-2)) {
        digits=abs(round(log10(qvalue[5])))+1
        value.labels=round(tick.at,digits=digits)
      } else {
        value.labels=round(tick.at,digits=2)
      }
    } else {
      if (is.null(tick.at)) {
        cat("tick.at should be specified if value.labels are defined.","\n")
      }
    }
    axis(side=2,lwd=0.8,cex.axis=1,at=tick.at,labels=value.labels)
    abline(h=ylim[1])
    mtext(side=3, title, line=1, cex=title.cex)
    if (!is.null(ylab)) {
      mtext(side=2,ylab,line=3,cex=2)
    }

    if (!is.null(indlist)) {
      if (text) {
        if (is.null(textwhat)) {
          text(x[indlist], value[indlist], indlist, col = colmat[indlist],
               cex = indcex)
        } else {
          value[indlist] = seq(min(value[indlist]), max(value[indlist]),
                               length.out = length(indlist))
          text(x[indlist], value[indlist], textwhat, col = colmat[indlist],
               cex = indcex)
        }
      } else {
        points(x[indlist],value[indlist],col=colmat[indlist],pch=19,cex=indcex)
      }
    }
    segments(x0=min(x),x1=max(x),y0=qvalue[3],y1=qvalue[3],col=box.col,lwd=5)
    segments(x0=min(x),x1=max(x),y0=qvalue[2],y1=qvalue[2],col=box.col,lwd=2)
    segments(x0=min(x),x1=max(x),y0=qvalue[4],y1=qvalue[4],col=box.col,lwd=2)
    segments(x0=min(x),x1=min(x),y0=qvalue[2],y1=qvalue[4],col=box.col,lwd=2)
    segments(x0=max(x),x1=max(x),y0=qvalue[2],y1=qvalue[4],col=box.col,lwd=2)

    if (print.outliers) {
      seg.range=max(x)-min(x)
      if (fence[1]>=ylim[1]) {
        segments(x0=min(x)+0.2*seg.range,x1=max(x)-0.2*seg.range,
                 y0=fence[1],y1=fence[1],col=box.col,lwd=2)
      }
      segments(x0=0.5,x1=0.5,y0=qvalue[2],y1=max(ylim[1],fence[1]),col=box.col,lwd=2)
      if (fence[2]<=ylim[2]) {
        segments(x0=min(x)+0.2*seg.range,x1=max(x)-0.2*seg.range,
                 y0=fence[2],y1=fence[2],col=box.col,lwd=2)
      }
      segments(x0=0.5,x1=0.5,y0=qvalue[4],y1=min(ylim[2],fence[2]),col=box.col,lwd=2)
    }
  }
  if (print.plot) {
    return(list(outliers=box.outliers,
                outliers.above=box.outliers.above,
                outliers.below=box.outliers.below,
                iqr=iqr,fence=fence))
  } else {
    return(list(outliers=box.outliers,
                outliers.above=box.outliers.above,
                outliers.below=box.outliers.below,
                iqr=iqr,fence=fence))
  }
}

