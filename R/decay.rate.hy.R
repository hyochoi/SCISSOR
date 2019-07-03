#'
#' @export
decay.rate.hy=function(Data,dist2end=NULL,robust=FALSE){
  # Calculate slopes & decay rates for every data point
  # Data = after applying center_data
  # dist2end = distance from the 3' end
  #            if null, use full transcript
  if (!is.null(dist2end)) {
    d = dim(data)[1];
    dist2end = min(d,dist2end);
    Data = Data[c((d-dist2end+1):d),];
    rm(d)
  }
  d = dim(Data)[1]; n = dim(Data)[2];
  lmcoef = matrix(0,ncol=2,nrow=n);
  x = 1:d ;
  for (case in 1:n){
    y = Data[,case];
    lmres = lm(y~x);
    lmcoef[case,] = coef(lmres);
  }
  rate = lmcoef[,2]*d
  return(list(slope=lmcoef[,2],rate=rate));
}
