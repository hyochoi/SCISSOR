#'
#' @export
lm.pval = function(x,y,true.sd=NULL,intercept=FALSE) {
  
  if (intercept) {
    mean.x = mean(x); mean.y = mean(y)
    x = x-mean.x; y = y-mean.y
    
    ssx = sum(x^2); lmdf = length(x)-2
    if (ssx>0) {
      b1 = sum(x*y)/ssx; b0 = mean.y-b1*mean.x
      lmcoef = c(b0,b1)
      if (is.null(true.sd)) {
        sd.y = sqrt(sum(((y+mean.y-b0-b1*x+b1*mean.x)^2)/lmdf))
      } else {
        sd.y = true.sd
      }
      
      b1.std = sd.y/sqrt(ssx)
      b0.std = sd.y*sqrt((1/length(x))+((mean.x^2)/ssx))
      lmstd = c(b0.std,b1.std)
    } else {
      b1 = 0; b0 = 0;
      lmcoef = c(b0,b1)
      if (is.null(true.sd)) {
        sd.y = sqrt(sum(((y+mean.y-b0-b1*x+b1*mean.x)^2)/lmdf))
      } else {
        sd.y = true.sd
      }
      lmstd = c(0,0)
    }
  } else {
    lmdf = length(x)-1
    ssx = sum(x^2)
    if (ssx>0) {
      b1 = sum(x*y)/ssx; b0 = 0
      lmcoef = c(b0,b1)
      if (is.null(true.sd)) {
        sd.y = sd(y-b1*x)
      } else {
        sd.y = true.sd
      }
      b1.std = sd.y/sqrt(ssx)
      b0.std = 1
      lmstd = c(b0.std,b1.std)
    } else {
      b1 = 0; b0 = 0;
      lmcoef = c(b0,b1)
      if (is.null(true.sd)) {
        sd.y = sqrt(sum(((y+mean.y-b0-b1*x+b1*mean.x)^2)/lmdf))
      } else {
        sd.y = true.sd
      }
      lmstd = c(0,0)
    }
  }
  
  if (sd.y>0) {
    tval = lmcoef/lmstd
    pval = 2*pt(q=abs(tval),df=lmdf,lower.tail=FALSE)
  } else {
    tval = c(0,0); pval = c(1,1);
  }
  return(list(lmcoef=lmcoef,lmstd=lmstd,tval=tval,pval=pval))
}
