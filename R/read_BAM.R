#' Read BAM files into R
#'
#' Get coverage pileup from BAM files for the provided regions.
#'
#' @param BAMfiles a list of BAM file names
#' @param caseIDs a list of case IDs in order of BAMfiles
#' @param gaf
#' @param symbol a gene symbol (gene name)
#' @param regions genomic regions to read formatted as "chr1:1-100,200-300:+". See also \code{\link{build_gaf}}.
#' @param outputType type of intronic region that will be included in output,
#'   with choices "whole_intron", "part_intron", or "only_exon"; the second is
#'   the default.
#' @param strand.specific indicator to get only "+" strand pileup
#' @examples
#' regions="chr1:1-100,200-300:+"
#' start_time = Sys.time()
#' countPileup = read_BAM(BAMfiles=BAMfiles,caseIDs=caseIDs,regions=regions)
#' end_time = Sys.time()
#' end_time - start_time
#'
#' @import Rsamtools GenomicRanges
#' @export
read_BAM = function(BAMfiles,caseIDs=NULL,gaf=NULL,symbol=NULL,regions=NULL,outputType="part_intron",strand.specific=FALSE,...) {

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

  covPileup = sapply(BAMfiles,function(x) read_aBAM(BAM=x,regions=new.regions,strand.specific=strand.specific,...))
  rownames(covPileup) = allPos
  colnames(covPileup) = caseIDs
  if (strnd=="+") {
    return(covPileup)
  } else {
    return(covPileup[rev(1:nrow(covPileup)),])
  }
}

#' @import Rsamtools GenomicRanges
#'
read_aBAM = function(BAM,regions=NULL,strand.specific=FALSE,...) {
  bf = Rsamtools::BamFile(BAM)
  ## prepare GRanges
  chr = strsplit(regions,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(regions,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(regions,":")[[1]][3]
  df = GenomicRanges::GRanges(chr,IRanges(start=as.numeric(strtend[,1]), end=as.numeric(strtend[,2])),strnd)

  # https://www.rdocumentation.org/packages/Rsamtools/versions/1.24.0/topics/pileup
  s_param = Rsamtools::ScanBamParam(which=df, what=c("pos"))
  # p_param = PileupParam(max_depth=1000, min_base_quality = 0, min_mapq = 0, min_nucleotide_depth = 0, min_minor_allele_depth = 0, distinguish_strands = FALSE, distinguish_nucleotides = FALSE, ignore_query_Ns = FALSE, include_deletions = FALSE, include_insertions = FALSE, left_bins = NULL, query_bins = NULL,  cycle_bins = NULL)
  if (!strand.specific) {
    res = pileup(bf, scanBamParam=s_param, pileupParam=PileupParam(distinguish_strands=F,distinguish_nucleotides=F,
                                                                   include_deletions = FALSE, include_insertions = FALSE,
                                                                   left_bins = NULL, query_bins = NULL,  cycle_bins = NULL,
                                                                   ...))
  } else {
    res = pileup(bf, scanBamParam=s_param, pileupParam=PileupParam(distinguish_strands=T,distinguish_nucleotides=F,
                                                                   include_deletions = FALSE, include_insertions = FALSE,
                                                                   left_bins = NULL, query_bins = NULL,  cycle_bins = NULL,
                                                                   ...))
    res = res[which(res$strand=="+"),]
  }
  ## pileup does not report zero coverage. manually put zero coverage to the output.
  ## will be deprecated or replaced to sth
  allPos = unlist(sapply(1:length(df), function(x) data.frame(df)[x,2]:data.frame(df)[x,3]))
  tmpDepth = rbind(res[,c("pos","count")],cbind(pos=allPos[!allPos %in% res$pos],count=rep(0,length(allPos[!allPos %in% res$pos]))))
  outpileup=tmpDepth[order(tmpDepth$pos),]
  rownames(outpileup)=1:dim(outpileup)[1]
  return(outpileup[,2])
}
