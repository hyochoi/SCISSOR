#' Build gene annotation
#'
#' Get gene annotation ranges (exons) for all genes from a GTF file.
#'
#' @param GTF.file GTF file name (with full directory). If NULL, UCSC annotation
#'   databases will be used.
#' @param output.file an output file name where the gene annotation ranges will be saved
#'
#' @examples
#' GTF.file = "genes.gtf"
#' regions = build_gaf(GTF.file=GTF.file)
#' write.table(x=regions,file="genes_gaf.txt",sep="\t",quote=F,row.names=F)
#'
#' @import GenomicRanges ballgown dplyr
#' @export
build_gaf <- function(GTF.file=NULL, output.file=NULL) {
  require(ballgown)
  require(GenomicRanges)
  require(dplyr)

  if (is.null(GTF.file)) {
    stop("GTF.file must be specified.")
  }
  gtfdf <- gffRead(GTF.file)
  cols <- c("seqname","start","end","strand")
  if (! "seqname" %in% colnames(gtfdf)) {
    if ("seqid" %in% colnames(gtfdf)) {
      cols <- c("seqid","start","end","strand")
    } else {
      stop("seqid is not available")
    }
  }
  ## Obtain for all genes
  get_exons <- function(gene_id,gtfdf,gtfGenes,cols) {
    tmp.gtfdf <- gtfdf[which(gtfGenes$gene_id==gene_id),]
    tmp.bed <- tmp.gtfdf[tmp.gtfdf$feature=="exon",cols]

    ## collect exons
    exondata <- data.frame(reduce(GRanges(tmp.bed)))
    exons <- paste(exondata$seqnames[1],paste(paste(exondata$start,exondata$end,sep="-"),
                                              collapse=","),
                   exondata$strand[1],sep=":")
    return(exons)
  }

  genes <- getAttributeField(gtfdf$attributes,"gene_name")
  genes <- sapply(genes,function(x) substr(x, 2, nchar(x)-1))
  geneids <- getAttributeField(gtfdf$attributes,"gene_id")
  geneids <- sapply(geneids,function(x) substr(x, 2, nchar(x)-1))
  gtfGenes <- data.frame(gene_name=genes, gene_id=geneids)

  ## Parallel approach
  gtfGenes_distinct = distinct(gtfGenes)
  cat(paste("\t Number of genes =",dim(gtfGenes_distinct)[1]),"\n")

  seqsplit <- rep(1:(round(dim(gtfGenes_distinct)[1]/1000)+1),each=1000)
  gtfGenes_distinct <- data.frame(gtfGenes_distinct, split=seqsplit[1:dim(gtfGenes_distinct)[1]])
  gtfGenes_distinct_split <- split(x=gtfGenes_distinct,f=gtfGenes_distinct$split)
  exons_subset <- vector("list",length(gtfGenes_distinct_split))
  for (isp in 1:length(gtfGenes_distinct_split)) {
    gtfGenes_distinct_subset <- gtfGenes_distinct_split[[isp]]
    gtfGenes_subset <- gtfGenes[which(gtfGenes$gene_id %in% gtfGenes_distinct_subset$gene_id),]
    gtfdf_subset <- gtfdf[which(gtfGenes$gene_id %in% gtfGenes_distinct_subset$gene_id),]
    exons_subset[[isp]] <- data.frame(gtfGenes_distinct_subset[,c(1,2)],
                                      regions=sapply(gtfGenes_distinct_subset$gene_id,
                                                     function(x) get_exons(gene_id=x,gtfdf=gtfdf_subset,
                                                                           gtfGenes=gtfGenes_subset,cols=cols)))
    cat(paste(c(paste(rep("=",isp),collapse = ""),
                paste(c(100*round(isp/length(gtfGenes_distinct_split),digits=2),"%")))),"\n")
  }
  all.regions <- do.call("rbind",exons_subset)
  if (! is.null(output.file)) {
    write.table(x=all.regions,
                file=output.file,
                sep="\t",quote=F,row.names=F)
    cat(paste(output.file,"has been created"),"\n")
  } else {
    write.table(x=all.regions,
                file=paste0(dirname(GTF.file),"/SCISSOR_gaf.txt"),
                sep="\t",quote=F,row.names=F)
    cat(paste(paste0(dirname(GTF.file),"/SCISSOR_gaf.txt"),"has been created"),"\n")
  }
  return(all.regions)
}
