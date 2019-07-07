#' This function is a wrapper of the necessary functions to detect RNA-seq shape
#' changes at a single gene in a single command.
#'
#' @param pileup a pile-up data matrix obtained from bam files, i.e. base-level
#'   raw counts. Columns are samples and rows are genomic positions.
#' @param exon an exon annotation for the given gene.
#' @param siglev a significance level for detecting outliers. Default is 1e-4.
#' @param logshiftVal a pseudo count added to raw counts before the logarithmic
#'   transformation. Default is NULL and automatically selects the pseudo count
#'   based on the implemented algorithm.
#' @param windowSize a window size used in the local shape change detection.
#' @param inputType type of intronic region contained in pileup, with choices
#'   "whole_intron", "part_intron", or "only_exon"; the second is the default.
#' @param outputType type of intronic region that will be included in output,
#'   with choices "whole_intron", "part_intron", or "only_exon"; the second is
#'   the default.
#' @param plotNormalization logical, whether to plot variation with respect to
#'   overall expression with fitted curve. Default is FALSE
#' @param reducedReturn
#'
#' @export
Scissor = function(pileup, exon, gene=NULL,
                   siglev=1e-4, logshiftVal=NULL, windowSize=50,
                   inputType="part_intron", outputType="part_intron",
                   plotNormalization=TRUE,
                   reducedReturn=FALSE) {

  ##---------------------------------------
  ##          Pre-process pileup
  ##---------------------------------------
  data.process = process_data(pileup=pileup,exon=exon,
                              logshiftVal=logshiftVal,
                              inputType=inputType,
                              outputType=outputType,
                              plotNormalization=plotNormalization,gene=gene)

  logshiftVal = data.process$logshiftVal  # log shift parameter
  exonset = data.process$dai$epm  # exon boundaries in data
  cat(paste0("     Log shift parameter used  = ",logshiftVal),"\n")

  ##---------------------------------------
  ##          Detect outliers
  ##---------------------------------------
  ## Step 1: Glocal shape change detection
  GSCout = miscGlobal(inputData=data.process$normalizedData,
                      siglev=siglev,ADcutoff=3,
                      PCnum=NULL,maxPCnum=10,
                      reducedReturn=reducedReturn)

  cat(paste0("     # Global Shape Changes identified  = ",length(GSCout$SC)),"\n")

  ## Step 2: Local shape change detection
  LSCout = miscLocal(miscGlobalResult=GSCout,
                       countData=data.process$countData,exonset=exonset,
                       siglev=siglev,cutoff=NULL,ADcutoff=3,
                       windowSize=windowSize,reducedReturn=reducedReturn)
  cat(paste0("     # Local Shape Changes identified   = ",length(LSCout$SC)),"\n")

  outputObject = list(SC=c(GSCout$SC,LSCout$SC), globalSC=GSCout$SC, localSC=LSCout$SC,
                      GSCout=GSCout, LSCout=LSCout,
                      countData=data.process$countData,
                      logData=data.process$logData,
                      normalizedData=data.process$normalizedData,
                      dai=data.process$dai,msf=data.process$msf,logshiftVal=logshiftVal)
  class(outputObject) = append(class(outputObject),"ScissorOutput")
  return(outputObject)
}
