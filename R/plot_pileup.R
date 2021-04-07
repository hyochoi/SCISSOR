#' Plot pileup profiles
#'
#' \code{plot_pileup} plots read counts of pileup coverage at base resolution.
#'
#' @param Pileup a pileup object with rows for genomic positions and columns for
#'   samples
#' @param Ranges a list of information for genomic regions.
#'   \code{\link{get_Ranges}}.
#' @param cases a vector specifying which samples are to be plotted. Sample IDs
#'   matching with the caseIDs in \code{\link{read_BAM}} or
#'   \code{\link{read_pileup}} can be used. Also, the column numbers can be
#'   simply used. e.g. \code{cases}=c(1,4,100) will plot the 1st, 4th, 100th
#'   columns (i.e. samples) in the \code{Pileup}. If NULL, all samples with
#'   plotted. Default is NULL.
#' @param logcount a numeric value indicating the number used in the logarithmic
#'   transformation. When your data is the log-transformed one, this parameter
#'   will print the raw read counts in y-axis labels instead of the log-scale.
#'   For example, if your data is \code{Pileup} = log10(read depth+c)-log10(c),
#'   "c" is the log count. Default is NULL, and y-axis will print labels based
#'   on the scale of the \code{Pileup}. If your data is not log-scale, use NULL.
#' @param plot.meanpileup a logical value indicating whether the mean coverage
#'   should be drawn. Default is TRUE.
#' @param col.pileup colors for the pileup samples
#' @param col.meanpileup a color for the mean coverage. Defult is "darkgrey"
#' @param print.ranges a logical value indicating whether genomic ranges and
#'   locus ranges (\code{\link{get_Ranges}}) should be printed in x-axis.
#' @param ... other graphical parameters passed to \code{points}.
#'
#' @keywords
#'
#' @examples
#' \donotrun{
#' plot_pileup(Pileup=ScissorOutput$logData,Ranges=Ranges,cases=2,logcount=10)
#' }
#'
#' @import RColorBrewer
#' @export
plot_pileup = function(Pileup,Ranges,cases=NULL,logcount=NULL,
                       plot.meanpileup=TRUE,
                       col.pileup=NULL,col.meanpileup="grey",
                       main=NULL,cex.main=1.2,
                       print.ranges=TRUE,
                       xlim=NULL,ylim=NULL,xlab=NULL,ylab=NULL,...) {

  ##  % needed variables
  if (missing(Pileup)) {
    stop("Pileup is missing")
  }
  if (missing(Ranges)) {
    stop("Genomic ranges should be needed (See get_Ranges)")
  }
  samplist = colnames(Pileup)
  n = ncol(Pileup)
  exons = matrix(Ranges$lRanges[,c(2,3)],ncol=2)
  if (is.null(cases)) cases = 1:n
  if (!is.numeric(cases)) {
    caseIDs=cases
    cases=which(samplist %in% caseIDs)
    if (length(cases)==0) {
      stop("No sample matches with the given case IDs.")
    }
  }

  require("RColorBrewer")
  candicol1 = c(brewer.pal(9,"Pastel1")[6], # candidate colors for exonic regions
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
  candicol2 = brewer.pal(12,"Set3") # candidiate colors for regions with shape changes
  candicol3 = brewer.pal(8,"Pastel2") # candidiate colors for regions with shape changes
  candicol = c(candicol2,candicol3);
  exon.col = candicol1[9]

  # Set plot parameters
  if (is.null(xlim)) xlim = c(0,nrow(Pileup))
  if (is.null(ylim)) ylim = yaxis.hy(Pileup)
  if (is.null(xlab)) xlab = "Genomic positions"
  if (is.null(ylab)) ylab = "Read depth"
  if (is.null(main)) {
    if (length(cases)>1) {
      main = paste0(Ranges$Gene)
    } else if (length(cases)==1) {
      main = paste0(Ranges$Gene," | sample #",cases," (ID:",samplist[cases],")")
    }
  }
  # pileup colors
  if (is.null(col.pileup)) {
    if (length(cases)>10) {
      x.mda = svd(Pileup) ;
      projmat = diag(x.mda$d)%*%t(x.mda$v) ;
      projmat[1,] = -projmat[1,] ;
      # colset = rainbow(n=length(projmat[1,]), start=0, end=0.756)[rank(-projmat[1,])]
      colset = colorRampPalette(brewer.pal(10, "Spectral"))(n)[rank(-projmat[1,])]
    } else {
      colset = rep("black",n)
    }
  } else {
    if (length(col.pileup)==n) {
      colset=col.pileup
    } else {
      colset = rep("black",n)
      colset[cases] = rep(col.pileup,length=length(cases))
    }
  }

  # Start plotting
  if (print.ranges) {
    par(mar=c(3.2,3.5,3,2))
  } else {
    par(mar=c(3,3.5,3,2))
  }
  meanPileup = apply(Pileup, 1, median) ;
  plot(meanPileup, type='l', lty=2, lwd=0.5, ylim=c(min(0,ylim[1]),ylim[2]),
       xlim=xlim, axes=F, ylab=NA, xlab=NA, xaxs="i",yaxs="i", col="white") ;
  for (i in 1:nrow(exons)){
    polygon(x=c(rep(exons[i,1],2),rep(exons[i,2],2)),y=c(-10000,(max(Pileup)+10000),(max(Pileup)+10000),-10000),col=exon.col,border=NA) ;
  }
  abline(v=exons[,1],lty=1,col="lightyellow3",lwd=0.1) ;
  abline(v=exons[,2],lty=1,col="lightyellow3",lwd=0.1) ;
  title(main, cex.main=cex.main,font.main=1,line=0.5);

  if (print.ranges) {
    x.tick.at = c(1,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
    x.labels.l = c(Ranges$Gene,Ranges$lRanges[2:nrow(Ranges$lRanges),2],max(Ranges$lRanges))
    x.labels.c = c(Ranges$Gene,Ranges$cRanges[2:nrow(Ranges$lRanges),1],max(Ranges$cRanges))
    x.labels.g = c(Ranges$chr,Ranges$gRanges[2:nrow(Ranges$gRanges),2],max(Ranges$gRanges))

    exon.tick.at = c(1,apply(Ranges$lRanges[,c(2,3)],1,mean))
    exon.labels = c(Ranges$Gene,paste("E",1:dim(Ranges$lRanges)[1],sep=""))
    axis(side=1, tck=-0.01, at=x.tick.at, labels=NA, col.ticks="darkgrey") ;
    # axis(side=1, lwd=0, line=-1, cex.axis=0.8,col.axis="darkgrey",
    #      at=x.tick.at,labels=x.labels.c) ;
    axis(side=1, lwd=0, line=-1, cex.axis=0.8,col.axis="darkgrey",
         at=exon.tick.at,labels=exon.labels) ;
    axis(side=1, lwd=0, line=-0.1, cex.axis=0.8,col.axis="darkgrey",
         at=x.tick.at,labels=x.labels.g) ;
    mtext(side=1, xlab, line=2, cex=1) ;
  } else {
    mtext(side=1, xlab, line=1, cex=1)
  }
  if (!is.null(logcount)) {
    if (logcount==1) {
      labels = c(1,5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    } else {
      labels = c(5,10,50,100,300,500,1000,2000,5000,10000,15000,20000,30000)
    }
    tick.at = log10(labels+logcount)-log10(logcount);
  } else {
    tick.at = NULL;
    labels = TRUE;
  }
  axis(side=2, tck=-0.02, at=tick.at, col.ticks="darkgrey",las=1,
       labels=labels,lwd=0,line=-0.8,cex.axis=0.8,col.axis="darkgrey")
  mtext(side=2, ylab, line=2, cex=1) ;

  box(lwd=1.5)
  if (plot.meanpileup){
    points(meanPileup, type='l', lty=1, lwd=2, col=col.meanpileup) ;
  }
  for (case in cases) {
    points(Pileup[,case], type='l', lty=1, col=colset[case], ...) ;
  }
}
