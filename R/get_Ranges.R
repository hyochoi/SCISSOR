#' Genomic ranges and locus ranges
#'
#' Get genomic ranges and locus ranges from the regions provided
#'
#' @param regions genomic regions formatted as chr1:1-100,200-300:+"
#' @param outputType type of intronic region that will be included in output,
#'   with choices "whole_intron", "part_intron", or "only_exon"; the second is
#'   the default.
#' @examples
#' regions="chr17:7571720-7573008,7573927-7574033,7576525-7576657,7576853-7576926,7577019-7577155,7577499-7577608,7578177-7578289,7578371-7578554,7579312-7579590,7579700-7579721,7579839-7579940:-"
#' Ranges=get_Ranges(regions=regions)
#'
#' @import BiocManager Rsamtools
#' @export
get_Ranges = function(Gene=NULL,regions,outputType="part_intron") {
  require(Rsamtools)

  if (missing(regions)) {
    stop("regions must be specified")
  }
  if (is.null(Gene)) {
    warning("The gene name is missing")
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
  intron.len = ceiling(len.intron.hy(exon=regions)*0.5);
  ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
  lRanges0=ep.new$ep
  gRanges0=cbind(strtend.num[,1]-(lRanges0[,2]-lRanges0[,1]),strtend.num,strtend.num[,2]+(lRanges0[,4]-lRanges0[,3]))
  new.regions=paste0(chr,":",paste(apply(gRanges0[,c(1,4)],1,function(x) paste(x,collapse="-")),collapse=","),":",strnd)

  if (strnd=="+") {
    gRanges=gRanges0
    lRanges=lRanges0
    colnames(gRanges)=c("ip.start","e.start","e.end","ip.end") # e:exon; ip: intronic part;
    rownames(gRanges)=paste("exon",1:nrow(gRanges),sep="")
    colnames(lRanges)=c("ip.start","e.start","e.end","ip.end")
    rownames(lRanges)=paste("exon",1:nrow(gRanges),sep="")
  } else {
    gRanges=gRanges0[rev(1:nrow(gRanges0)),rev(1:4)]
    lRanges=max(lRanges0)-lRanges0[rev(1:nrow(lRanges0)),rev(1:4)]+1
    colnames(gRanges)=c("ip.start","e.start","e.end","ip.end") # e:exon; ip: intronic part;
    rownames(gRanges)=paste("exon",1:nrow(gRanges),sep="")
    colnames(lRanges)=c("ip.start","e.start","e.end","ip.end")
    rownames(lRanges)=paste("exon",1:nrow(gRanges),sep="")
  }
  rm(gRanges0,lRanges0)
  output=list(Gene=Gene,gRanges=gRanges,lRanges=lRanges,chr=chr,strand=strnd,regions=regions,new.regions=new.regions)
  output
}
