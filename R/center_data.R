#' Centering data
#'
#' This function finds a robust mean vector along the gene length using a trimmed mean. Also, it provides
#' a new data object
#'
#' @export
center_data = function(inputData, ...) {
  msfres = scale.factor(X=inputData, ...)
  data.center = msfres$mean.vec
  msf=msfres$msf
  msf2=rep(1,length(msf))
  msf2[which(msf>1)]=1/msf[which(msf>1)]
  outdata = sweep(msfres$NewX,2,msf2,"*")
  return(list(outputData=outdata,msf=msf,data.center=data.center))
}
