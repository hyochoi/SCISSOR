#'
#' @export
pca2misc = function(scoremat,eigenvalue,dimension,
                    siglev=1e-4,ADcutoff=3,
                    PCnum=NULL,maxPCnum=10,
                    reducedReturn=FALSE) {
  # Version 3: Add choosing PC directions
  # Computes the Stahel-Donoho outlyingness of every element in x
  # x can be univariate or multivariate
  # Assumes that dim(x) = n x d if multivariate
  #
  n = ncol(scoremat); d = dimension;
  eigenval = eigenvalue;
  projmat = scoremat;

  if (length(which(eigenval>0))<10) {
    K=0
  } else {
    # Choose the number of spikes
    if (is.null(PCnum)){
      result=PCnum.gamma.hy(eigenval=eigenval,d=d,n=n,eps=(1/d));
      K = max(1,min(result$K,maxPCnum));
    } else {
      K = PCnum;
    }
  }
  # Get OS statistic
  if (K==0) {
    # cat(paste("!!!!!!!     ","No directions included","    !!!!!!!"),"\n");
    OS = rep(0,n);
    OSpval = rep(1,n);
    NPS = NULL;  # Projection Score
    if (reducedReturn) {
      NPS = NULL;
    }
    return(list(OS=OS,OSpval=OSpval,NPS=NPS,
                outliers=NULL,outOS=NULL,outliers.sort=NULL,outOS.sort=NULL,
                K=K));
  } else {
    x = projmat[1:K,];   # projmat with the chosen directions

    if (is.null(dim(x))){
      ADstat = ADstatWins.hy(x);
      if (ADstat > ADcutoff) {
        OS = rep(0,n);
        OSpval = rep(1,n);
        NPS = NULL;  # Projection Score
        directions = NULL;
        return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                    outliers=NULL,outOS=NULL,outliers.sort=NULL,outOS.sort=NULL,
                    K=K));
      } else {
        OS = abs(x-median(x))/mad(x);
        OSpval = 1-pchisq(q=OS^2,df=K);
        NPS = matrix(rep(OS,n),ncol=n);
        directions = matrix(rep(x,n),ncol=n);
        if (reducedReturn) {
          NPS = directions = NULL;
        }
        out = which(OSpval<siglev); outOS = OS[out];
        out.order = order(outOS,decreasing=TRUE);
        out.sort = out[out.order]; outOS.sort = outOS[out.order];
        return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                    outliers=out,outOS=outOS,outliers.sort=out.sort,outOS.sort=outOS.sort,
                    K=K));
      }# if (ADstat > ADcutoff)
    } else {
      x = t(x);
      ndir = 300*K;
      A = generdir(x,ndir=ndir) # generates `ndir' directions (ndir by K)

      # Add more potential directions
      subprojmat = t(x);
      Scov = (subprojmat%*%t(subprojmat))/n;
      B0 = solve(Scov)%*%subprojmat;   # K by n
      B0 = B0[,-which(apply(B0,2,FUN=function(x){sqrt(sum(x^2))})<1e-10)];
      B1 = sweep(B0,2,apply(B0,2,FUN=function(x){sqrt(sum(x^2))}),"/")  # normalized K by n
      A = rbind(A,t(B1))  # ndir by K

      Y = x %*% t(A) # project x onto A (n by ndir)
      ADstat = apply(Y,2,ADstatWins.hy);
      indir = which(ADstat<ADcutoff);

      if (length(indir)>1) {
        A = A[indir,]; Y = Y[,indir];
        out_temp = apply(X=Y,MARGIN=2,
                         FUN= function(t) (t-median(t))/mad(t)); # n by length(indir)
        indexmax = apply(abs(out_temp),1,which.max);
        NPS = out_temp[,indexmax];
        directions=t(A[indexmax,]);

        neg.index = which(diag(NPS)<0);
        if (length(neg.index)>0) {
          NPS[,neg.index] = -NPS[,neg.index];
          directions[,neg.index] = -directions[,neg.index];
        }
        OS=diag(NPS)
        OSpval = 1-pchisq(q=OS^2,df=K);
      } else {
        OS = rep(0,n);
        OSpval = rep(1,n);
        NPS = directions= NULL;  # Projection Score
      }

      if (reducedReturn) {
        NPS = directions = NULL;
      }
      out = which(OSpval<siglev); outOS = OS[out];
      out.order = order(outOS,decreasing=TRUE);
      out.sort = out[out.order]; outOS.sort = outOS[out.order];
      return(list(OS=OS,OSpval=OSpval,NPS=NPS,directions=directions,
                  outliers=out,outOS=outOS,outliers.sort=out.sort,outOS.sort=outOS.sort,
                  K=K));
    }
  }
}
