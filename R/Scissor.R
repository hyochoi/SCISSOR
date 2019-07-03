#' Performs the shape change detection in a single command
#'
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
Scissor = function(pileup, exon, siglev=1e-4, logshiftVal=NULL, windowSize=50,
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
                              plotNormalization=plotNormalization)

  logshiftVal = data.process$logshiftVal  # log shift parameter
  exonset = dai$epm  # exon boundaries in data
  cat(paste0("     Log shift parameter used  = ",logshiftVal),"\n")

  ##---------------------------------------
  ##          Detect outliers
  ##---------------------------------------
  ## Step 1: Glocal shape change detection
  GSCout = miscGlobal(inputData=data.process$normalizedData,siglev=siglev,
                        PCnum=NULL,eps=NULL,maxPCnum=10,
                        ADcutoff=3,reducedReturn=reducedReturn)

  globalSC = GSCout$outliers.sort
  cat(paste0("     # Global Shape Changes identified  = ",length(globalSC)),"\n")

  ## Step 2: Local shape change detection
  LSCout = miscLocal(miscGlobalobj=GSCout,
                       countData=data.process$countData,exonset=exonset,
                       siglev=siglev,cutoff=NULL,ADcutoff=3,
                       windowSize=windowSize,reducedReturn=reducedReturn)
  localSC = LSCout$outliers.sort
  cat(paste0("     # Local Shape Changes identified   = ",length(localSC)),"\n")

  return(list(SC=c(globalSC,localSC), globalSC=globalSC, localSC=localSC,
              GSCout=GSCout, LSCout=LSCout,
              countData=data.process$countData,
              logData=data.process$logData,
              normalizedData=data.process$normalizedData,
              dai=data.process$dai,msf=data.process$msf,logshiftVal=logshiftVal))
}
