#'
#' @param eps tuning parameter for estimating the number of PCs. If NULL
#'   (default), eps = (1/nrow(data))
#'
#' @import nloptr
#' @export
PCnum.gamma.hy = function(eigenval,d,n,eps=1e-4){
  # Third version - consider jth eigenvalue as a spike only if it is larger than a criteria
  #                   using (j+2)th to last eigenvalues.
  #                   Previous version used (j+1)th to last eigenvalues
  # Determine the number of spikes when the underlying distribution for PSD is a truncated gamma(shape,rate)
  require(nloptr);
  lss.gamma = function(theta,m,xi=(1-1e-4)){
    # m = ESD moments in vector
    k = theta[1]; scl = 1/theta[2]; # k=shape; lam=scale
    u = qgamma(xi,shape=k,scale=scl)
    sq.sum = (m[1]*xi - (k*scl)*pgamma(u,shape=(k+1),scale=scl))^2 + (m[2]*xi - (scl^2)*k*(k+1)*pgamma(u,shape=(k+2),scale=scl))^2
    return(sq.sum)
  }

  stiel <- function(u,d,n,eigvals){
    #	Sample Stieltjes transform
    inner_sum <- 0 ;
    m <- length(eigvals)
    for (i in 1:m){
      inner_sum <- inner_sum + 1/(eigvals[i]-u) ;
    }
    if (d < n){
      fval <- inner_sum/n - (1-(d/n))/u ;
    } else {
      fval <- inner_sum/n ;
    }
    return(fval) ;
  }

  phi.gamma = function(x,y,shape=1,rate=1,lower=0,upper=2){
    in.integrate = function(t){
      (t/(x-t))*dgamma(x=t,shape=shape,rate=rate)
    }
    int = integrate(f=in.integrate,lower=lower,upper=upper)$value
    multi.scale=pgamma(upper,shape=shape,rate=rate)
    return(x + y*x*int*(1/multi.scale))
  }
  ##############################
  d = d-1;
  u.thrsh = c(); theta.hat <- matrix(,ncol=2,nrow=0); iter = c();
  jc = 1 ; jp = 0;  # jc: spike to be tested at current step; jp: spike tested at previous step
  while (jc > jp) {
    iter = c(iter,jc);
    jp = jc;
    a = c(); b = c(); # a: ESD moments; b: PSD moments
    d = d-1; y = d/n;
    a[1] = sum(eigenval[(jc+2):length(eigenval)])/(d);
    a[2] = sum(eigenval[(jc+2):length(eigenval)]^2)/(d);
    b[1] = a[1]; b[2] = a[2] - y*(a[1]^2);
    theta.0 = c();   # theta[1]=shape, theta[2]=rate
    theta.0[2] = a[1]/(a[2]-(1+y)*(a[1]^2)); theta.0[1] = a[1]*theta.0[2];

    obj.f = function(x){ lss.gamma(theta=x,m=b,xi=(1-eps))}
    S1 = bobyqa(theta.0,obj.f,lower=c(0,0),upper=c(1e+5,1e+5))
    theta.1 = S1$par

    upper = qgamma((1-eps),shape=theta.1[1],rate=theta.1[2])
    sn = function(u){stiel(u=u,d=d,n=n,eigvals=eigenval[(jc+2):n])}
    alpha = -1/sn(eigenval[jc]);

    expr = function(x){phi.gamma(x=x,y=y,shape=theta.1[1],rate=theta.1[2],upper=upper)}
    thrsh.alpha = optimize(expr,interval=c(upper,1000))$minimum  # threshold for alpha (true eigenvalue)
    if (expr(thrsh.alpha) > eigenval[jc]){
      jc = jc;
    } else {
      jc = jc + 1;
      u.thrsh = c(u.thrsh,expr(thrsh.alpha));
      theta.hat = rbind(theta.hat,theta.1);
    }
  }
  if (jc > 1) {
    return(list(K=(jc-1),theta=theta.hat[(jc-1),],theta.his=theta.hat[1:(jc-1),],threshold=u.thrsh,iter=iter));
  } else {
    return(list(K=0,theta=theta.1,theta.his=theta.1,threshold=expr(thrsh.alpha),iter=iter))
  }
}
