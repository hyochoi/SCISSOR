#' Read BAM files into R
#'
#' Get coverage pileup from BAM files for the provided regions.
#'
#' @param BAMfiles
#' @param caseIDs
#' @param gaf
#' @param symbol
#' @param regions genomic regions formatted as chr1:1-100,200-300:+"
#' @param outputType type of intronic region that will be included in output,
#'   with choices "whole_intron", "part_intron", or "only_exon"; the second is
#'   the default.
#' @examples
#' manifest=read.table(file="~/Desktop/HNSC/base/HNSC_RNA_manifest_datastore.txt")
#' manifest=manifest[1:10,]
#' caseIDs=as.character(manifest[,1])
#' BAMfiles=apply(manifest,1,function(x) paste(c("~/Desktop/HNSC_BAM",strsplit(x=as.character(x[2]),split="/")[[1]][c(6,7,8)]),collapse="/"))
#' regions="chr17:7571720-7573008,7573927-7574033,7576525-7576657,7576853-7576926,7577019-7577155,7577499-7577608,7578177-7578289,7578371-7578554,7579312-7579590,7579700-7579721,7579839-7579940:-"
#' start_time = Sys.time()
#' countPileup = read_BAM(BAMfiles=BAMfiles,caseIDs=caseIDs,regions=regions)
#' end_time = Sys.time()
#' end_time - start_time
#'
#' @import BiocManager Rsamtools
#' @export
read_BAM = function(BAMfiles,caseIDs=NULL,gaf=NULL,symbol=NULL,regions=NULL,outputType="part_intron",...) {

  require(Rsamtools)

  if (missing(BAMfiles)) {
    stop("BAM file path is missing")
  }
  if (is.null(caseIDs)) {
    caseIDs = 1:length(BAMfiles)
  }
  if (is.null(regions)) {
    if (is.null(symbol)) {
      stop("Either symbol or regions should be specified")
    } else {
      if (is.null(gaf)) {
        stop("gaf path is missing")
      } else {
        regions = NULL # should be fixed.
      }
    }
  }

  Ranges=get_Ranges(regions=regions,outputType=outputType)
  new.regions=Ranges$new.regions
  chr = strsplit(new.regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(new.regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(new.regions,":")[[1]][3]
  strtend.num=matrix(as.numeric(strtend),ncol=2)
  allPos = unlist(sapply(1:nrow(strtend.num), function(x) strtend.num[x,1]:strtend.num[x,2]))

  covPileup = sapply(BAMfiles,function(x) read_aBAM(BAM=x,regions=new.regions))
  rownames(covPileup) = allPos
  colnames(covPileup) = caseIDs
  if (strnd=="+") {
    return(covPileup)
  } else {
    return(covPileup[rev(1:nrow(covPileup)),])
  }
}

read_aBAM = function(BAM,regions=NULL,...) {
  require(Rsamtools)
  bf = BamFile(BAM)
  ## prepare GRanges
  chr = strsplit(regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(regions,":")[[1]][3]
  df = GRanges(chr,IRanges(start=as.numeric(strtend[,1]), end=as.numeric(strtend[,2])),strnd)

  s_param = ScanBamParam(which=df, what=c("pos"))
  # p_param = PileupParam(max_depth=1000, min_base_quality = 0, min_mapq = 0, min_nucleotide_depth = 0, min_minor_allele_depth = 0, distinguish_strands = FALSE, distinguish_nucleotides = FALSE, ignore_query_Ns = FALSE, include_deletions = FALSE, include_insertions = FALSE, left_bins = NULL, query_bins = NULL,  cycle_bins = NULL)
  res = pileup(bf, scanBamParam=s_param, pileupParam=PileupParam(distinguish_strands=F,distinguish_nucleotides=F,...))
  ## pileup does not report zero coverage. manually put zero coverage to the output.
  ## will be deprecated or replaced to sth
  allPos = unlist(sapply(1:length(df), function(x) data.frame(df)[x,2]:data.frame(df)[x,3]))
  tmpDepth = rbind(res[,c("pos","count")],cbind(pos=allPos[!allPos %in% res$pos],count=0))
  outpileup=tmpDepth[order(tmpDepth$pos),]
  rownames(outpileup)=1:dim(outpileup)[1]
  return(outpileup[,2])
}
