#' Genomic ranges and locus ranges
#'
#' Get genomic ranges and locus ranges from the regions provided
#'
#' @param regions genomic regions formatted as chr1:1-100,200-300:+"
#' @examples
#' regions="chr17:7571720-7573008,7573927-7574033,7576525-7576657,7576853-7576926,7577019-7577155,7577499-7577608,7578177-7578289,7578371-7578554,7579312-7579590,7579700-7579721,7579839-7579940:-"
#' Ranges=get_Ranges(regions=regions)
#'
#' @import BiocManager Rsamtools
#' @export
get_Ranges = function(regions) {
  require(Rsamtools)
  chr = strsplit(regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(regions,":")[[1]][3]
  df = GRanges(chr,IRanges(start=as.numeric(strtend[,1]), end=as.numeric(strtend[,2])),strnd)
  strtend.num=matrix(as.numeric(strtend),ncol=2)

  intron.len = ceiling(len.intron.hy(exon=regions)*0.5);
  ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
  lRanges=ep.new$ep
  gRanges=cbind(strtend.num[,1]-(lRanges[,2]-lRanges[,1]),strtend.num,strtend.num[,2]+(lRanges[,4]-lRanges[,3]))
  new.Ranges=paste0(chr,":",paste(apply(gRanges[,c(1,4)],1,function(x) paste(x,collapse="-")),collapse=","),":",strnd)

  if (strnd=="+") {
    colnames(gRanges)=c("ip.start","e.start","e.end","ip.end") # e:exon; ip: intronic part;
    rownames(gRanges)=paste("exon",1:nrow(gRanges),sep="")
    colnames(lRanges)=c("ip.start","e.start","e.end","ip.end")
    rownames(lRanges)=paste("exon",1:nrow(gRanges),sep="")
  } else {
    colnames(gRanges)=rev(c("ip.start","e.start","e.end","ip.end"))
    rownames(gRanges)=rev(paste("exon",1:nrow(gRanges),sep=""))
    colnames(lRanges)=rev(c("ip.start","e.start","e.end","ip.end") )
    rownames(lRanges)=rev(paste("exon",1:nrow(gRanges),sep=""))
  }

  output=list(gRanges=gRanges,lRanges=lRanges,chr=chr,strand=strnd,regions=regions,new.Ranges=new.Ranges)
  output
}
