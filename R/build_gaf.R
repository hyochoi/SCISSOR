#' Build gene annotation
#'
#' Get genome ranges (exons) for a given gene from UCSC or a given GTF file.
#'
#' @param Gene character(1) identifying the gene symbol, e.g. \code{"TP53"}
#' @param GTF.file GTF file name (with full directory). If NULL, UCSC annotation
#'   databases will be used.
#' @param hg.ref either "hg19" or "hg38", which will be used for UCSC annotation
#'   databases. By default, hg19 is used.
#'
#' @examples
#' build_gaf(Gene="TP53")
#'
#' @import BiocManager GenomicRanges refGenome
#' @export
build_gaf = function(Gene,GTF.file=NULL,hg.ref=c("hg19","hg38")) {
  if (missing(Gene)) {
    stop("Gene symbol is missing")
  }
  if (length(Gene)>1) {
    warning("More than one gene was provided. The first gene will be used.")
    Gene=Gene[1]
  }

  if (is.null(GTF.file)) {
    exons=c()
  #   require(org.Hs.eg.db)
  #   hg.ref = match.arg(hg.ref, choices=c("hg19","hg38"))
  #
  #   if (hg.ref=="hg19") {
  #     require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  #     txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
  #   } else {
  #     require(TxDb.Hsapiens.UCSC.hg38.knownGene)
  #     txdb = TxDb.Hsapiens.UCSC.hg38.knownGene
  #   }
  #
  #   geneid = select(org.Hs.eg.db, keys=Gene, columns=c("ENTREZID"),keytype="SYMBOL")
  #   ebg = exonsBy(txdb, by="gene") ## Once loaded, can be used at multiple instances
  #   exondata = data.frame(reduce(ebg[which(names(ebg)==as.character(geneid[2]))]))
  #   exons = paste(exondata$seqnames[1],paste(paste(exondata$start,exondata$end,sep="-"),collapse=","),
  #                 exondata$strand[1],sep=":")
  } else {
    require(refGenome)
    require(GenomicRanges)

    ens = ensemblGenome()
    # GTF
    read.gtf(ens,GTF.file, useBasedir=FALSE) ## read.gtf does not accept gtf.gz: gunzip it.
    tmp = getGtf(extractByGeneName(ens,Gene))
    # tmp = getGtf(extractByGeneId(ens,geneid$GENEID))
    tmp.bed = tmp[tmp$feature=="exon",c("seqid","start","end","strand")]
    exondata = data.frame(reduce(GRanges(tmp.bed)))
    if (sum(!grepl(pattern="chr",exondata$seqnames))>0) {
      exondata$seqnames[(!grepl(pattern="chr",exondata$seqnames))] = paste("chr",exondata$seqnames[(!grepl(pattern="chr",exondata$seqnames))],sep="")
    }
    exons = paste(exondata$seqnames[1],paste(paste(exondata$start,exondata$end,sep="-"),collapse=","),
                  exondata$strand[1],sep=":")
  }

  return(exons)
}

