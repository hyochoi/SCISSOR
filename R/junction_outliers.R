#'
#' @export
get_SpliceRatio = function(Ranges,JSR.table,JSR.matrix) {

  lRanges = Ranges$lRanges
  nexons = nrow(lRanges)
  n = ncol(JSR.matrix)
  junctions.g = matrix(as.numeric(split_junction(as.character(JSR.table$junctions.g))),ncol=2,byrow=T)
  junctions.l = matrix(as.numeric(split_junction(as.character(JSR.table$junctions.l))),ncol=2,byrow=T)
  junction.start = unique(junctions.g[,1])

  ## Exon splice ratio
  exon.splice.ratio=matrix(0,ncol=ncol(JSR.matrix),nrow=nrow(JSR.matrix))
  for (i in 1:length(junction.start)) {
    index.temp=which(junctions.g[,1]==junction.start[i])
    junction.sum = apply(matrix(JSR.matrix[index.temp,],ncol=n),2,sum);

    if (max(junction.sum)>0) {
      for (j in 1:n) {
        if (junction.sum[j]>10) {
          exon.splice.ratio[index.temp,j] = JSR.matrix[index.temp,j]/junction.sum[j]
        }
      }
      for (j in 1:n) {
        if (junction.sum[j]<=10) {
          for (k in index.temp) {
            exon.splice.ratio[k,j] = median(exon.splice.ratio[k,])
          }
        }
      }
    }
  }
  colnames(exon.splice.ratio) = colnames(JSR.matrix)

  ## Intron splice ratio
  intron.splice.ratio=matrix(0,ncol=ncol(JSR.matrix),nrow=2*(nexons-1))
  for (ie in 1:(nexons-1)) {
    ## 5' intron splice ratio
    k = 2*(ie-1) + 1
    coverage.depth0=rawData[(lRanges[ie,3]+10),]
    index.temp=which(junctions.l[,1]==lRanges[ie,3])
    junction.sum = apply(matrix(JSR.matrix[index.temp,],ncol=n),2,sum)
    coverage.depth = coverage.depth0 + junction.sum

    intron.splice.ratio[k,which(coverage.depth>10)]=1-(junction.sum[which(coverage.depth>10)]/coverage.depth[which(coverage.depth>10)])
    intron.splice.ratio[k,which(coverage.depth<=10)]=median(intron.splice.ratio[k,which(coverage.depth>10)])
    # kdeplot.hy(intron.splice.ratio[k,],xlim=c(0,1),main=paste(ie,"| 5'"))

    ## 3' intron splice ratio
    k = 2*(ie-1) + 2
    coverage.depth0=rawData[(lRanges[(ie+1),2]-10),]
    index.temp=which(junctions.l[,2]==lRanges[(ie+1),2])
    junction.sum = apply(matrix(JSR.matrix[index.temp,],ncol=n),2,sum)
    coverage.depth = coverage.depth0 + junction.sum

    intron.splice.ratio[k,which(coverage.depth>10)]=1-(junction.sum[which(coverage.depth>10)]/coverage.depth[which(coverage.depth>10)])
    intron.splice.ratio[k,which(coverage.depth<=10)]=median(intron.splice.ratio[k,which(coverage.depth>10)])
    # kdeplot.hy(intron.splice.ratio[k,],xlim=c(0,1),main=paste(ie,"| 3'"))
  }
  rownames(intron.splice.ratio) = rep(c(5,3),times=(nexons-1))
  colnames(intron.splice.ratio) = colnames(JSR.matrix)

  return(list(exon.splice.ratio=exon.splice.ratio,intron.splice.ratio=intron.splice.ratio))
}

#'
#' @export
find_intron_outliers_ie = function(ie,Ranges,JSR.matrix,JSR.table,splice.ratio,
                                   siglev=1e-04,cutoff.a=0.1) {
  ## input = i (which row), JSR.matrix, JSR.table, splice.ratio
  ## Set boxplot cutoff
  q4=qnorm(0.75)
  iqr.range=((qnorm(1-(siglev/2))/q4)-1)/2 # cutoff for choosing outliers
  lRanges = Ranges$lRanges
  nexons = nrow(lRanges)

  junctions.g = matrix(as.numeric(split_junction(as.character(JSR.table$junctions.g))),ncol=2,byrow=T)
  junctions.l = matrix(as.numeric(split_junction(as.character(JSR.table$junctions.l))),ncol=2,byrow=T)

  k5 = 2*(ie-1) + 1
  k3 = 2*(ie-1) + 2
  ## 5' intron retention
  coverage.depth0=rawData[(lRanges[ie,3]+10),]
  index.temp=which(junctions.l[,1]==lRanges[ie,3])
  i = index.temp[which(junctions.l[index.temp,2]==lRanges[(ie+1),2])]
  junction.sum = apply(matrix(JSR.matrix[index.temp,],ncol=n),2,sum)
  coverage.depth = coverage.depth0 + junction.sum

  a0=boxplot.hy(splice.ratio[k5,],robust=T,print.plot=F,
                ylim=c(0,max(splice.ratio)),
                title=junctions.g[ie],title.cex=0.9,
                indlist=8,indlist.col="red",
                outliers.range=iqr.range,outliers.col="grey")
  # kdeplot.hy(splice.ratio[ie,],indlist=a0$outliers.above,xlim=c(0,1),text=T,main=ie)
  outliers5.tmp = outliers5.L1 = outliers5.L2 = c()
  if ((length(which(splice.ratio[k5,]>cutoff.a))<0.1*n) & (length(a0$outliers.above)>0)) {
    outliers5.tmp = a0$outliers.above[which(splice.ratio[k5,a0$outliers.above]>=cutoff.a)]
    if (length(outliers5.tmp)>0) {
      outliers5.L1 = outliers5.tmp[which((splice.ratio[k5,outliers5.tmp]-(1-ES_cutoff_fn(coverage.depth[outliers5.tmp])))>0)] # level 1 outliers
      outliers5.L2 = outliers5.tmp[which(! outliers5.tmp %in% outliers5.L1)]
    }
  }

  ## 3' intron retention
  coverage.depth0=rawData[(lRanges[(ie+1),2]-10),]
  index.temp=which(junctions.l[,2]==lRanges[(ie+1),2])
  junction.sum = apply(matrix(JSR.matrix[index.temp,],ncol=n),2,sum)
  coverage.depth = coverage.depth0 + junction.sum

  a0=boxplot.hy(splice.ratio[k3,],robust=T,print.plot=F,
                ylim=c(0,max(splice.ratio)),
                title=junctions.g[ie],title.cex=0.9,
                indlist=8,indlist.col="red",
                outliers.range=iqr.range,outliers.col="grey")
  # kdeplot.hy(splice.ratio[ie,],indlist=a0$outliers.above,xlim=c(0,1),text=T,main=ie)
  outliers3.tmp = outliers3.L1 = outliers3.L2 = c()
  if ((length(which(splice.ratio[k3,]>cutoff.a))<0.1*n) & (length(a0$outliers.above)>0)) {
    outliers3.tmp = a0$outliers.above[which(splice.ratio[k3,a0$outliers.above]>=cutoff.a)]
    if (length(outliers3.tmp)>0) {
      outliers3.L1 = outliers3.tmp[which((splice.ratio[k3,outliers3.tmp]-(1-ES_cutoff_fn(coverage.depth[outliers3.tmp])))>0)] # level 1 outliers
      outliers3.L2 = outliers3.tmp[which(! outliers3.tmp %in% outliers3.L1)]
    }
  }
  outliers5 = unique(c(outliers5.L1,outliers5.L2))
  outliers3.L1 = outliers3.L1[which(! outliers3.L1 %in% outliers5)]
  outliers3.L2 = outliers3.L2[which(! outliers3.L2 %in% outliers5)]

  output = rbind(matrix(cbind(outliers5.L1,
                              rep("intron",times=length(outliers5.L1)),
                              rep("5",times=length(outliers5.L1)),
                              rep("Level_1",times=length(outliers5.L1)),
                              round(splice.ratio[k5,outliers5.L1],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers5.L1)),
                              as.matrix(JSR.table[rep(i, each = length(outliers5.L1)), ])),ncol=(6+ncol(JSR.table))),
                 matrix(cbind(outliers5.L2,
                              rep("intron",times=length(outliers5.L2)),
                              rep("5",times=length(outliers5.L2)),
                              rep("Level_2",times=length(outliers5.L2)),
                              round(splice.ratio[k5,outliers5.L2],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers5.L2)),
                              as.matrix(JSR.table[rep(i, each = length(outliers5.L2)), ])),ncol=(6+ncol(JSR.table))),
                 matrix(cbind(outliers3.L1,
                              rep("intron",times=length(outliers3.L1)),
                              rep("3",times=length(outliers3.L1)),
                              rep("Level_1",times=length(outliers3.L1)),
                              round(splice.ratio[k3,outliers3.L1],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers3.L1)),
                              as.matrix(JSR.table[rep(i, each = length(outliers3.L1)), ])),ncol=(6+ncol(JSR.table))),
                 matrix(cbind(outliers3.L2,
                              rep("intron",times=length(outliers3.L2)),
                              rep("3",times=length(outliers3.L2)),
                              rep("Level_2",times=length(outliers3.L2)),
                              round(splice.ratio[k3,outliers3.L2],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers3.L2)),
                              as.matrix(JSR.table[rep(i, each = length(outliers3.L2)), ])),ncol=(6+ncol(JSR.table))))
  if (dim(output)[1]>0) {
    colnames(output) = c("Outlier","Region","Sign","Level","VJF","Junction.name",colnames(JSR.table))
    output = data.frame(output)
    output$JV.class = rep("IR",dim(output)[1])
    output$Region.tag = paste("I",sapply(split_junction(as.character(output$LBE.position))[1,],FUN=function(t){strsplit(strsplit(t,":")[[1]][1],"exon")[[1]][2]}),sep="")
  } else {
    output = NULL
  }
  return(output)
}

#'
#' @export
find_exon_outliers_i = function(i,JSR.matrix,JSR.table,splice.ratio,siglev=1e-04,
                                cutoff.a=0.1,cutoff.b=0.9) {
  ## input = i (which row), JSR.matrix, JSR.table, splice.ratio
  ## Set boxplot cutoff
  q4=qnorm(0.75)
  iqr.range=((qnorm(1-(siglev/2))/q4)-1)/2 # cutoff for choosing outliers

  junctions.g = matrix(as.numeric(split_junction(as.character(JSR.table$junctions.g))),ncol=2,byrow=T)
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
  outliers.tmp1 = outliers.tmp2 = outliers.tmp1.L1 = outliers.tmp1.L2 = outliers.tmp2.L1 = outliers.tmp2.L2 = c()
  # if ((length(which(splice.ratio[i,]<cutoff.b))<0.1*n) & (length(a0$outliers.below)>0)) {
  #   ## This is the case of possible exon skipping / cryptic events / SV
  #   outliers.tmp1 = a0$outliers.below[which(splice.ratio[i,a0$outliers.below]<cutoff.b)]
  #   outliers.tmp1.L1 = outliers.tmp1[which((splice.ratio[i,outliers.tmp1]-ES_cutoff_fn(junction.sum[outliers.tmp1]))<0)]
  #   # outliers.tmp1 = outliers.tmp1[sapply(outliers.tmp1,FUN=function(k) {(JSR.matrix[i,k]<(0.8*median(JSR.matrix[which(JSR.matrix[,k]>5),k])))})]
  #   outliers.tmp1.L2 = outliers.tmp1[which(! outliers.tmp1 %in% outliers.tmp1.L1)]
  # }
  if ((length(which(splice.ratio[i,]>cutoff.a))<0.1*n) & (length(a0$outliers.above)>0)) {
    ## This is the case of possible alternative exon / cryptic events / SV
    outliers.tmp2 = a0$outliers.above[which(splice.ratio[i,a0$outliers.above]>=cutoff.a)]
    outliers.tmp2.L1 = outliers.tmp2[which(sapply(outliers.tmp2,FUN=function(k) { (JSR.matrix[i,k]>max(median(JSR.matrix[which(JSR.matrix[,k]>=5),k])*0.2,10))})==T)]
    outliers.tmp2.L2 = outliers.tmp2[which(! outliers.tmp2 %in% outliers.tmp2.L1)]
  }

  output = rbind(matrix(cbind(outliers.tmp1.L1,
                              rep("exon",times=length(outliers.tmp1.L1)),
                              rep("-",times=length(outliers.tmp1.L1)),
                              rep("Level_1",times=length(outliers.tmp1.L1)),
                              round(splice.ratio[i,outliers.tmp1.L1],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers.tmp1.L1)),
                              as.matrix(JSR.table[rep(i, each = length(outliers.tmp1.L1)), ])),
                        ncol=(6+ncol(JSR.table))),
                 matrix(cbind(outliers.tmp1.L2,
                              rep("exon",times=length(outliers.tmp1.L2)),
                              rep("-",times=length(outliers.tmp1.L2)),
                              rep("Level_2",times=length(outliers.tmp1.L2)),
                              round(splice.ratio[i,outliers.tmp1.L2],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers.tmp1.L2)),
                              as.matrix(JSR.table[rep(i, each = length(outliers.tmp1.L2)), ])),
                        ncol=(6+ncol(JSR.table))),
                 matrix(cbind(outliers.tmp2.L1,
                              rep("exon",times=length(outliers.tmp2.L1)),
                              rep("+",times=length(outliers.tmp2.L1)),
                              rep("Level_1",times=length(outliers.tmp2.L1)),
                              round(splice.ratio[i,outliers.tmp2.L1],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers.tmp2.L1)),
                              as.matrix(JSR.table[rep(i, each = length(outliers.tmp2.L1)), ])),
                        ncol=(6+ncol(JSR.table))),
                 matrix(cbind(outliers.tmp2.L2,
                              rep("exon",times=length(outliers.tmp2.L2)),
                              rep("+",times=length(outliers.tmp2.L2)),
                              rep("Level_2",times=length(outliers.tmp2.L2)),
                              round(splice.ratio[i,outliers.tmp2.L2],digits=3),
                              rep(rownames(JSR.table)[i],times=length(outliers.tmp2.L2)),
                              as.matrix(JSR.table[rep(i, each = length(outliers.tmp2.L2)), ])),
                        ncol=(6+ncol(JSR.table))))

  if (dim(output)[1]>0) {
    colnames(output) = c("Outlier","Region","Sign","Level","VJF","Junction.name",colnames(JSR.table))
    output = data.frame(output)
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

