#' Provide the winsorized Anderson-Darling test statistic
#'
#' This function winsorizes and calculates the A-D test statistic.
#' Based on the MATLAB function ADStatQF by Qing Feng.
#'
#' @param x an input vector
#'
#' @export
ADstatWins.hy=function(x){
  # Winsorize and calculate Anderson Darling test statistic
  n = length(x);
  mabd = mean(abs(x-median(x)));

  if (mabd==0){
    xfinal = rep(0,n);
  } else {
    xk = (x-median(x))/mabd ;
    # Winsorise data
    p95 = 1.0972;   # Extreme value inverse with p=95%, mu=0,sigma=1
    an = (2*log(n))^(-0.5);
    bn = sqrt(2*log(n)-log(log(n))-log(4*pi));
    L = p95*an+bn;
    xk[which(xk>L)] = rep(L,length(which(xk>L)));
    xk[which(xk<(-L))] = rep((-L),length(which(xk<(-L))));

    # Re-standardize
    xfinal = (xk-mean(xk))/sd(xk);
  }

  # Calculate evaluation stat
  ADval = ADstat.hy(xfinal);
  return(ADval);
}
