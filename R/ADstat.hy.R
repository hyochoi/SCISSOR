#' Provide the Anderson-Darling test statistic of standard normal distribution
#'
#' This function is based on ADStatQF function in MATLAB by Qing Feng.
#'
#' @param x an input vector
#'
#' @export
ADstat.hy=function(x){
  n = length(x);
  if (n < 7){
    return(print("Sample size must be greater than 7"));
  } else {
    xs = sort(x);
    f = pnorm(xs,mean(xs),sd(xs));
    i = 1:n ;
    S = sum(((2*i-1)/n)*(log(f)+log(1-rev(f))));
    ADval = -n-S;
    return(ADval);
  }
}
