#'
check_normals = function(normal.data,exon,out1=NULL,out2=NULL,
                         data.process=NULL,global.out=NULL,local.out=NULL,
                         input.type="part_intron",output.type="part_intron",intron.len=NULL,
                         plot.g=FALSE,plot.l=FALSE,colmatT=NULL,cutoff=NULL,main=NULL, ...) {

  if (is.null(out1)) {
    out1=global.out$outliers.sort
  }
  if (is.null(out2)) {
    out2=local.out$outliers
  }
  main0 = main; rm(main)
  dataS = data.process$dataS
  logshift.val = data.process$logshift.val
  mean.vec = data.process$data.center

  n = ncol(dataS); d = nrow(dataS); n2 = ncol(GeneDataN)
  # log transform
  normal.annotation=annotate_pileup(pileup=GeneDataN,exon=exon,input.type=input.type,
                                    output.type=output.type,intron.len=intron.len)

  datalog.N = log10(normal.annotation$out.pileup+logshift.val)-log10(logshift.val)
  msf.N = as.vector(t(mean.vec)%*%datalog.N/(sum(mean.vec^2)));
  dataS.N = datalog.N - mean.vec%*%t(msf.N)

  # Prepare datasets for two-step method
  K = max(c(global.out$PCsubset,global.out$rm.PCdir))
  dataS.Ng = global.out$pca$dirmat[,1:K]%*%(t(global.out$pca$dirmat[,1:K])%*%dataS.N) # low-rank matrix
  dataS.Nl = dataS.N - dataS.Ng # residual matrix

  dataS.Tg = dataS - global.out$resmat
  dataS.Tl = global.out$resmat
  # collect MOD
  MODn.nps.g = MODn.g = nps.g.max = nps.g.whichmax = out1n = out1n.which = NULL
  MODn.nps.l = MODn.l = nps.l.max = nps.l.whichmax = out2n = out2n.which = NULL
  if (length(out1)>0) {
    MODn.nps.g = matrix(0,ncol=length(out1),nrow=(n+n2))
    dataS.TN=cbind(dataS.Tg,dataS.Ng)
    for (j in 1:length(out1)) {
      currout=out1[j]
      MODn.nps.g[,j] = pd.rate.hy(as.vector(t(global.out$MOD[,currout])%*%dataS.TN))
      if (plot.g) {
        if (is.null(colmatT)) {
          colmat=rep("grey",(n+n2));
          colmat[c((n+1):(n+n2))] = "green"; colmat[out1] = "red"
        } else {
          colmat=c(colmatT,rep("green",n2))
        }
        if (is.null(main0)) {
          main=paste("Global -",out1[j])
        }
        if (is.null(cutoff)) {
          cutoff=global.out$cutoff
        }
        kdeplot.hy(MODn.nps.g[,j],main=main,colmat=colmat, ...)
        abline(v=cutoff,col="red")
      }
    }
    nps.g.max=apply(matrix(MODn.nps.g[c((n+1):(n+n2)),],ncol=length(out1)),1,max)
    nps.g.whichmax=apply(matrix(MODn.nps.g[c((n+1):(n+n2)),],ncol=length(out1)),1,which.max)
    out1n=which(nps.g.max>global.out$cutoff)
    out1n.which=out1[nps.g.whichmax[out1n]]
    MODn.g=global.out$MOD[,out1n.which]
    rm(dataS.TN)
  }
  if (length(out2)>0) {
    MODn.nps.l = matrix(0,ncol=length(out2),nrow=(n+n2))
    dataS.TN=cbind(dataS.Tl,dataS.Nl)
    for (j in 1:length(out2)) {
      currout=out2[j]
      MODn.nps.l[,j] = pd.rate.hy(as.vector(t(local.out$MOD[,currout])%*%dataS.TN))
      if (plot.l) {
        if (is.null(colmatT)) {
          colmat=rep("grey",(n+n2));
          colmat[c((n+1):(n+n2))] = "green"; colmat[out2] = "blue"
        } else {
          colmat=c(colmatT,rep("green",n2))
        }
        if (is.null(main0)) {
          main=paste("Local -",out2[j])
        }
        if (is.null(cutoff)) {
          cutoff=local.out$cutoff
        }
        kdeplot.hy(MODn.nps.l[,j],main=main,colmat=colmat, ...)
        abline(v=cutoff,col="red")
      }
    }
    nps.l.max=apply(matrix(MODn.nps.l[c((n+1):(n+n2)),],ncol=length(out2)),1,max)
    nps.l.whichmax=apply(matrix(MODn.nps.l[c((n+1):(n+n2)),],ncol=length(out2)),1,which.max)
    out2n = which(nps.l.max>local.out$cutoff)
    out2n.which=out2[nps.l.whichmax[out2n]]
    MODn.l=local.out$MOD[,out2n.which]
  }

  # output = out1n, out1n.which, out2n, out2n.which, MOD.nps.g, MOD.nps.l
  return(list(out1n=out1n,out1n.which=out1n.which,
              out2n=out2n,out2n.which=out2n.which,
              MODn.nps.g=MODn.nps.g,MODn.g=MODn.g,
              MODn.nps.l=MODn.nps.l,MODn.l=MODn.l,
              datalog.N=datalog.N,msf.N=msf.N,
              dataS.N=dataS.N,dataS.Ng=dataS.Ng,dataS.Nl=dataS.Nl))
}
