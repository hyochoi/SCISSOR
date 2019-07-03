#'
#' @export
SlogADstat.hy = function(X,shift.val=1){
  # log transform
  datalog = log10(X + shift.val) - log10(shift.val) ;
  # calculate the median scale factor
  msf = scale.factor(X=datalog,average="mean",trim=0.1,adjval=NULL)$msf

  ADval = ADstatWins.hy(msf);
  return(list(msf=msf,ADval=ADval));
}
