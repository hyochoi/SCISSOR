#' Processing data
#'
#' @export
process_data = function(pileup,exon,
                        inputType="whole_intron",outputType="part_intron",
                        logshiftVal=NULL,
                        plotNormalization=FALSE, ...) {

  # Annotate pileup data
  data.annotation = annotate_pileup(pileup=pileup,exon=exon,
                             inputType=inputType,outputType=outputType)
  dai = data.annotation$dai # data annotation information
  exonset = dai$epm
  intron.len = dai$intron.len

  # Log transform data
  data.log = logtransform_data(inputData=data.annotation$new.pileup,logshiftVal=logshiftVal)

  # Center and normalize data
  data.normalized = normalize_data(inputData=data.log$outputData,
                                   pileup=data.annotation$new.pileup,
                                   exonset=exonset,
                                   makePlot=plotNormalization, ...)

  readconstr=20
  exon.rm=c(); intron.rm=c();
  if (nrow(exonset)>1) {
    for (ie in 1:nrow(exonset)) {
      b.exon=c(exonset[ie,2]:exonset[ie,3])
      if (max(apply(data.annotation$new.pileup[b.exon,],1,max))<readconstr) {
        exon.rm=c(exon.rm,ie)
      }
    }

    for (ie in 1:(nrow(exonset)-1)) {
      b.exon1=c(exonset[ie,2]:exonset[ie,3])
      b.exon2=c(exonset[ie+1,2]:exonset[ie+1,3])
      q.exon=apply(cbind(apply(data.annotation$new.pileup[b.exon1,],2,FUN=function(x){quantile(x,probs=0.5)}),
                         apply(data.annotation$new.pileup[b.exon2,],2,FUN=function(x){quantile(x,probs=0.5)})),1,max)

      b.intron=c((exonset[ie,3]+6):(exonset[ie+1,2]-6)) # remove each boundary of length 4
      lcarea=length(b.intron)
      if (lcarea<22) { # remove introns less than 30 (22+4+4)
        intron.rm=c(intron.rm,ie)
      } else {
        q.intron=apply(cbind(apply(data.annotation$new.pileup[b.intron[6:ceiling(lcarea/3)],],2,FUN=function(x){quantile(x,probs=0.75)}),
                             apply(data.annotation$new.pileup[b.intron[ceiling(lcarea/3):(2*ceiling(lcarea/3))],],2,FUN=function(x){quantile(x,probs=0.75)}),
                             apply(data.annotation$new.pileup[b.intron[(2*ceiling(lcarea/3)):(lcarea-6)],],2,FUN=function(x){quantile(x,probs=0.75)})),1,max)

        if ((max(apply(data.annotation$new.pileup[b.intron,],1,max))<readconstr) &
            (length(which((q.intron>0.2*q.exon) & (q.exon>10)))==0)) {
          intron.rm=c(intron.rm,ie)
        }
      }
      # cat(ie,"\n")
    }

    new.exonset=exonset
    for (ie in 1:(nrow(exonset)-1)) {
      if (exonset[ie+1,2]==exonset[ie+1,1]) {
        baselen=floor(0.5*(exonset[ie+1,1]-exonset[ie,3]))
        new.exonset[ie,4]=exonset[ie,3]+baselen
        new.exonset[ie+1,1]=exonset[ie,3]+baselen+1
      }
    }
    # new.exonset

    if (length(exon.rm)>0) {
      for (ie in exon.rm) {
        data.normalized$outputData[c(new.exonset[ie,1]:new.exonset[ie,4]),]=0
      }
    }
    if (length(intron.rm)>0) {
      for (ie in intron.rm) {
        data.normalized$outputData[c((exonset[ie,3]+1):(exonset[ie+1,2]-1)),]=0
      }
    }
  }

  return(list(countData=data.annotation$new.pileup,
              logData=data.log$outputData,
              normalizedData=data.normalized$outputData,
              dai=dai,exonset=exonset,
              logshiftVal=data.log$logshiftVal,
              msf=data.normalized$msf,
              g1.offset=data.normalized$g1.offset,g2.offset=data.normalized$g2.offset,
              goodcase=data.normalized$goodcase,
              data.center=data.normalized$data.center))
}
