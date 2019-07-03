scale1StepM = function(x,precScale=1e-10) {
  # Computes the first step of an algorithm for
  # a scale M-estimator using the given rho function. 
  # The scatter is computed relative to zero.
  #
  x = x[!is.na(x)] # we always take out NAs
  n = length(x)
  if(n == 0) { return(0.0)
  } else {
    sigma0 = 1.4826*median(abs(x))
    if(sigma0 < precScale) { return(0.0)
    } else {
      rho = rhoHuber(x/sigma0)
      return(sigma0 * sqrt(sum(rho)*2/n))
    }
  }
}
