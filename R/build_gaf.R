#' Build gene annotation
#'
#' Get genome ranges (exons) for a given gene from a GTF file.
#'
#' @param Gene character(1) identifying the gene symbol for gene list, e.g. \code{"TP53"}
#' @param GTF.file GTF file name (with full directory). If NULL, UCSC annotation
#'   databases will be used.
#'
#' @examples
#' GTF.file = "genes.gtf"
#' regions = build_gaf(Gene=genelist[1:100],GTF.file=GTF.file)
#' write.table(x=regions,file="genes_gaf.txt",sep="\t",quote=F,col.names=F,row.names=F)
#'
#' @import GenomicRanges ballgown
#' @export
build_gaf = function(Gene,GTF.file=NULL) {
  require(ballgown)
  require(GenomicRanges)

  if (is.null(GTF.file)) {
    stop("GTF.file must be specified.")
  }
  gtfdf = gffRead(GTF.file)
  cols = c("seqname","start","end","strand")
  if (! "seqname" %in% colnames(gtfdf)) {
    if ("seqid" %in% colnames(gtfdf)) {
      cols = c("seqid","start","end","strand")
    } else {
      stop("seqid is not available")
    }
  }
  genes = getAttributeField(gtfdf$attributes,"gene_name")
  genes = sapply(genes,function(x) substr(x, 2, nchar(x)-1))

  ## Obtain for all genes
  get_exons = function(Gene,gtfdf,genes,cols) {
    tmp = gtfdf[which(genes==Gene),]
    tmp.bed = tmp[tmp$feature=="exon",cols]

    ## collect exons
    exondata = data.frame(reduce(GRanges(tmp.bed)))
    exons = paste(exondata$seqnames[1],paste(paste(exondata$start,exondata$end,sep="-"),collapse=","),
                  exondata$strand[1],sep=":")
    return(exons)
  }

  if (is.null(Gene)) {
    genelist = unique(genes)
    cat(paste("Gene is not specified and all genes in GTF.file will be processed."),"\n")
    cat(paste("\t Number of genes =",length(genelist)),"\n")
  } else {
    genelist = Gene
  }

  exons = sapply(genelist,function(x) get_exons(Gene=x,gtfdf=gtfdf,genes=genes,cols=cols))
  regions = data.frame(Genes=genelist,regions=exons)
  return(regions)
}
