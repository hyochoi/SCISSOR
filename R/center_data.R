#' Centering data
#'
#' This function finds a robust mean vector along the gene length using a trimmed mean. Also, it provides
#' a new data object
#'
#' @export
center_data = function(inputData, Ranges, average="mean",trim=0.1,adjval=NULL, ...) {

  d = dim(inputData)[1]
  nexons = nrow(Ranges$lRanges)
  mean.vec = apply(inputData,1,FUN=function(x){adj.center(x,average=average,trim=trim,adjval=adjval)})

  exonbaseMat = matrix(0,ncol=nexons,nrow=d)
  colnames(exonbaseMat) = rownames(Ranges$lRanges)
  for (ie in 1:nexons) {
    exonbaseMat[Ranges$lRanges[ie,2]:Ranges$lRanges[ie,3],ie] = mean.vec[Ranges$lRanges[ie,2]:Ranges$lRanges[ie,3]]/sum(mean.vec[Ranges$lRanges[ie,2]:Ranges$lRanges[ie,3]]^2)
  }
  msf.exons = t(exonbaseMat)%*%inputData
  msf = apply(msf.exons,2,median)

  NewX = inputData - mean.vec%*%t(msf)
  # msf2=rep(1,length(msf))
  # msf2[which(msf>1)]=1/msf[which(msf>1)]
  # NewX = sweep(NewX,2,msf2,"*")
  return(list(outputData=NewX,msf=msf,data.center=mean.vec))
}

center_data_v1 = function(inputData, robustLM=FALSE, ...) {
  msfres = scale.factor(X=inputData, robustLM=robustLM, ...)
  data.center = msfres$mean.vec
  msf=msfres$msf
  msf2=rep(1,length(msf))
  msf2[which(msf>1)]=1/msf[which(msf>1)]
  outdata = sweep(msfres$NewX,2,msf2,"*")
  return(list(outputData=outdata,msf=msf,data.center=data.center))
}
