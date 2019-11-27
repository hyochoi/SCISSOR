#'
#' @export
get_PO = function(X,siglev=1e-4,NormCutoff=3,canDir=NULL,numSearch=300) {

  if (is.null(X)) {
    stop("object is null")
  }

  if (is.null(dim(X))) { # X is a vector
    n = length(X)
    ADstat = ADstatWins.hy(X);
    if (ADstat < NormCutoff) {
      OS = abs(X-median(X))/mad(X)
      OSpval = 1-pchisq(q=OS^2,df=1)
      cutoff = qchisq(p=(1-siglev),df=1)
      NPS = matrix(rep(OS,n),ncol=n)
      directions = matrix(rep(X,n),ncol=n)
    } else {
      message("no direction satisfied the constraint")
      OS = rep(0,n)
      OSpval = rep(1,n)
      cutoff = 0
      NPS = directions = NULL
    }
  } else { # X is a matrix
    if (is.null(canDir)) {
      n = ncol(X); M = nrow(X);
      Scov = (X%*%t(X))/n
      d = length(which(eigen(Scov)$values>1e-13))
      x = t(X)
      ndir = numSearch*M
      A = generdir(x,ndir=ndir) # generates `ndir' directions (ndir by M)

      ## Add more potential directions
      if (M == d) {
        B0 = solve(Scov)%*%X   # M by n
        B0 = B0[,-which(apply(B0,2,FUN=function(x){sqrt(sum(x^2))})<1e-10)]
        B1 = sweep(B0,2,apply(B0,2,FUN=function(x){sqrt(sum(x^2))}),"/")  # normalized M by n
        A = rbind(A,t(B1))  # ndir by M
      }
      A = rbind(A,diag(M))
      ## Compute projection outlyingness
      Y = x %*% t(A) # project x onto A (n by ndir)
      ADstat = apply(Y,2,ADstatWins.hy) #
      indir = which(ADstat<NormCutoff)
      if (length(indir)>1) {
        A = A[indir,]; Y = Y[,indir]
        out_temp = apply(X=Y,MARGIN=2,
                         FUN= function(t) (t-median(t))/mad(t)) # n by length(indir)
        indexmax = apply(abs(out_temp),1,which.max)
        NPS = out_temp[,indexmax]
        directions=t(A[indexmax,])

        neg.index = which(diag(NPS)<0)
        if (length(neg.index)>0) {
          NPS[,neg.index] = -NPS[,neg.index];
          directions[,neg.index] = -directions[,neg.index];
        }
        OS=diag(NPS)
        OSpval = 1-pchisq(q=OS^2,df=d);
        cutoff = qchisq(p=(1-siglev),df=d);
      } else {
        message("no direction satisfied the constraint")
        OS = rep(0,n)
        OSpval = rep(1,n)
        cutoff = 0
        NPS = directions = NULL
      }
    }
  } # end of if (is.null(dim(X)))
  if (sum(OS)<1e-5) {
    out = NULL; out.sort = NULL
  } else {
    out = which(OSpval<siglev)
    out.sort = out[order(OS[out],decreasing=TRUE)]
  }
  return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,cutoff=cutoff,
              outliers=out,outliers.sort=out.sort))
}
