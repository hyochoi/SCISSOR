#' Genomic ranges and locus ranges
#'
#' This function returns genomic ranges and locus ranges in the needed format
#' for downstream SCISSOR analysis
#'
#' @param Gene character(1) identifying the gene symbol, e.g. \code{"TP53"}
#' @param regions character(1) describing genomic regions formatted as
#'   \code{"chr1:1-100,200-300:+"}
#' @param GTF.file GTF file name (with full directory). If NULL, UCSC annotation
#'   databases will be used. See \code{\link{build_gaf}.
#' @param hg.ref either "hg19" or "hg38", which will be used for UCSC annotation
#'   databases. By default, hg19 is used. See \code{\link{build_gaf}.
#' @param outputType character(1) of the type of intronic region that will be
#'   included in output, with choices "whole_intron", "part_intron", or
#'   "only_exon"; the second is the default.
#'
#' @return a \code{list} of key information for the regions provided: \itemize{
#'   \item Gene gene symbol \item gRanges ranges }
#'
#' @examples
#' regions="chr17:7571720-7573008,7573927-7574033,7576525-7576657,7576853-7576926,7577019-7577155,7577499-7577608,7578177-7578289,7578371-7578554,7579312-7579590,7579700-7579721,7579839-7579940:-"
#' Ranges=get_Ranges(Gene="TP53",regions=regions,outputType="part_intron")
#'
#' @import BiocManager Rsamtools
#' @export
get_Ranges = function(Gene=NULL,regions=NULL,GTF.file=NULL,hg.ref=c("hg19","hg38"),outputType="part_intron") {

  require(Rsamtools)
  if (!is.null(Gene)) {
    if (length(Gene)>1) {
      warning("More than one gene was provided. The first gene will be used.")
      Gene=Gene[1]
    }
  }
  if (is.null(regions)) {
    if (is.null(Gene)) {
      stop("Either Gene or regions must be specified")
    } else {
      hg.ref = match.arg(hg.ref, choices=c("hg19","hg38"))
      regions = build_gaf(Gene=Gene,GTF.file=GTF.file,hg.ref=hg.ref)
    }
  }
  chr = strsplit(regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(regions,":")[[1]][3]
  strtend.num=matrix(as.numeric(strtend),ncol=2)

  if (outputType=="whole_intron") {
    ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=NULL) ;
  } else if (outputType=="part_intron") {
    intron.len = ceiling(len.intron.hy(exon=regions)*0.5);
    ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
  } else if (outputType=="only_exon") {
    intron.len = ceiling(len.intron.hy(exon=regions)*0.5);
    ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=0) ;
  }
  lRanges0=ep.new$ep
  gRanges0=cbind(strtend.num[,1]-(lRanges0[,2]-lRanges0[,1]),strtend.num,strtend.num[,2]+(lRanges0[,4]-lRanges0[,3]))
  cRanges0=matrix(find.exon.hy(regions,is.intron=TRUE,num.intron=0)$ep[,c(2,3)],ncol=2)
  new.regions=paste0(chr,":",paste(apply(matrix(gRanges0[,c(1,4)],ncol=2),1,function(x) paste(x,collapse="-")),collapse=","),":",strnd)

  if (strnd=="+") {
    gRanges=gRanges0
    lRanges=lRanges0
    cRanges=cRanges0
  } else {
    gRanges=matrix(gRanges0[rev(1:nrow(gRanges0)),rev(1:4)],ncol=4)
    lRanges=matrix(max(lRanges0)-lRanges0[rev(1:nrow(lRanges0)),rev(1:4)]+1,ncol=4)
    cRanges=matrix(max(cRanges0)-cRanges0[rev(1:nrow(cRanges0)),2:1]+1,ncol=2)
  }
  rm(gRanges0,lRanges0,cRanges0)
  colnames(gRanges)=c("ip.start","e.start","e.end","ip.end") # e:exon; ip: intronic part;
  rownames(gRanges)=paste("exon",1:nrow(gRanges),sep="")
  colnames(lRanges)=c("ip.start","e.start","e.end","ip.end")
  rownames(lRanges)=paste("exon",1:nrow(gRanges),sep="")
  colnames(cRanges)=c("e.start","e.end")
  rownames(cRanges)=paste("exon",1:nrow(gRanges),sep="")

  output=list(Gene=Gene,gRanges=gRanges,lRanges=lRanges,cRanges=cRanges,chr=chr,strand=strnd,regions=regions,new.regions=new.regions)
  output
}
