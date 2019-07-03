#' Weighing data
#' @export
weigh_data = function(data,dai=NULL,method="mean",max.weight=2) {

  w = rep(1,nrow(data)); # initial weight
  outdata = data
  if (method=="equal") {
    outdata = outdata
  } else {
    epm = dai$epm; nexons = nrow(epm)
    epl = epm[,2]; epr = epm[,3]
    intron.len = dai$intron.len

    # Here, I gave a constraint for intron.len but should be more careful
    if (intron.len<10) {
      outdata = outdata
    } else {
      if (method=="max") {
        cutoff.length = 2*intron.len;
        max.intron.length = max(epl[2:nexons]-epr[1:(nexons-1)]-1) ;
        target.length = max.intron.length;
      } else if (method=="mean") {
        cutoff.length = floor(sum(epr-epl+1)/(nexons-1)) ;
        target.length = cutoff.length;
      }

      # Here, I gave a constraint for intron.len but should be more careful
      short.int = which(((epl[2:nexons]-epr[1:(nexons-1)]-1)>10) &
                          ((epl[2:nexons]-epr[1:(nexons-1)]-1)<cutoff.length));
      if (length(short.int)>0){
        #  Creat a weight matrix for short introns to have comparable variations.
        loc.w = c() ;
        for (kk in 1:length(short.int)){
          foo = (epm[short.int[kk],3]+1):(epm[(short.int[kk]+1),2]-1); # short intron base pairs
          if (length(foo)>10){
            loc.w = c(loc.w,foo) ;
            if (!max.weight) {
              w[foo] = sqrt(target.length/length(foo));
            } else {
              w[foo] = min(sqrt(target.length/length(foo)),2);
            }
          }
          #         print((2*num.intron[ii])/length(foo));
        }
        if (length(loc.w)>1){
          outdata[loc.w,] = diag(w[loc.w])%*%data[loc.w,];
        }
      }
    }
  }
  return(list(outdata=outdata,w=w))
}
