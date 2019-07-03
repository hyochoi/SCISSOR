#' Shift parameter
getShiftParam = function(X, param.grid=NULL, draw.plot=FALSE) {

  if (is.null(param.grid)){
    param.grid = seq(1,10,by=1)
  }
  res = apply(matrix(param.grid,ncol=1),1,
              FUN=function(x){SlogADstat.hy(X=X,shift.val=x)$ADval})
  optim.idx = which.min(res)
  optim.param = param.grid[optim.idx]
  optim.stat = res[optim.idx]

  if (is.infinite(optim.stat)) {
    warning("Minimum = Inf")
  }

  if (draw.plot){
    plot(param.grid,res,type="o",ylab="A-D statistic",xlab="Shift value")
    abline(v=optim.param,col="red");
    text(x=optim.param+(0.2*(max(param.grid)-min(param.grid))),
         y=(min(res)+0.8*(max(res)-min(res))),labels=round(optim.param,digits=3),col="red")
  }
  return(list(optim.param=optim.param,optim.stat=optim.stat))
}
