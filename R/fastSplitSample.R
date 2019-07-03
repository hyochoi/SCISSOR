fastSplitSample = function(x){
  # Centers sample by median, and divides in 2 equal halves.
  # Assumes that NAs have already been removed.
  # This function has time complexity O(n).
  #
  med = median(x) # takes only O(n) time
  x = x - med # centering
  n = length(x)  
  h = n %/% 2   # = integer part of n/2
  xa = x[x > 0] # length(xa) <= h
  xb = x[x < 0] # length(xa) <= h 
  xa = c(rep(0,(n - h - length(xa))),xa)
  xb = c(rep(0,(n - h - length(xb))),abs(xb)) # abs() ! 
  return(list(xa=xa,xb=xb,med=med))
}
