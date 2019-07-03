#' Detect local shape outliers based on a residual matrix
#'
#' This function identifies local shape variants by mining locally on/off abnormalities.
#'
#' @param resmat a residual matrix from PCA, a data matrix subtracted by a low-rank matrix.
#' @param rawmat a raw coverage matrix, dataI from \link{process_data}
#' @param exonset exon/intron annotation matrix from \link{annotate_pileup}
#' @param winlength the window size. Default is 100.
#' @param readconstr the minimum reads count required to be considered for on/off local shape changes.
#' Default is 10.
#'
#' @export
detect_onoffout = function(resmat,rawmat,exonset,
                           winlength=100,readconstr=10,
                           siglev=NULL,cutoff=NULL) {

  res_on = get_onstat(resmat=resmat,rawmat=rawmat,exonset=exonset,
                      readconstr=readconstr)
  res_off = get_offstat(resmat=resmat,rawmat=rawmat,exonset=exonset,
                        winlength=winlength,readconstr=readconstr)

  onoff_stat = apply(rbind(res_on$on_stat,res_off$off_stat),2,max)
  fivenumber = fivenum(onoff_stat)
  MOD = matrix(0,nrow=nrow(resmat),ncol=ncol(resmat))

  if (fivenumber[3]<=0) {
    cutoff=NULL; onoff_site = NULL;
    nout = 0;
  } else {
    onoff_i = apply(rbind(res_on$on_stat,res_off$off_stat),2,which.max)
    # 1 = on; 2 = off;
    whichon = which(onoff_i==1); whichoff = which(onoff_i==2);
    onoff_site = matrix(0,nrow=ncol(resmat),ncol=2)
    onoff_site[whichon,] = res_on$on_site[whichon,]
    onoff_site[whichoff,] = res_off$off_site[whichoff,]

    for (j in 1:ncol(resmat)) {
      temp_vec = rep(0,nrow(resmat))
      temp_site = onoff_site[j,]
      temp_vec[temp_site[1]:temp_site[2]] = 1/sqrt(temp_site[2]-temp_site[1]+1)
      if (onoff_i[j]==2) {
        temp_vec  = -temp_vec
      }
      MOD[,j] = temp_vec
    }

    if (is.null(cutoff)) {
      rsc=compScales(onoff_stat)
      cutoff = qnorm(p=(1-siglev))*rsc$sa + fivenumber[3]
    }
    nout = length(which(onoff_stat>cutoff))
  }

  if (nout>0) {
    onoffout = order(onoff_stat,decreasing=T)[1:nout]
    onout = order(res_on$on_stat,decreasing=T)[1:length(which(res_on$on_stat>cutoff))];
    offout = order(res_off$off_stat,decreasing=T)[1:length(which(res_off$off_stat>cutoff))];
  } else {
    onoffout = onout = offout = NULL
  }
  return(list(outliers=onoffout,on_outliers=onout,off_outliers=offout,MOD=MOD,
              cutoff=cutoff,
              onoff_stat=onoff_stat,onoff_site=onoff_site,
              res_on=res_on,res_off=res_off))
}
