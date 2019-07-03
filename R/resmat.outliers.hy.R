resmat.outliers.hy = function(resmat,GeneData,ep.new,seqdir="+",iqr.cutoff=5){
  ##  Outlier detection - 2nd step
  ##  Detect local outliers from residual matrix
  ##  Updated on September 16, 2017.
  ##  Hyo Young Choi
  
  cep.new = ep.new
  # if (seqdir=="-") {
  #   seqlen = ep.new$ep[nrow(ep.new$ep),4]
  #   cep.new[,1] = seqlen - rev(ep.new$ep[,4]) + 1;
  #   cep.new[,2] = seqlen - rev(ep.new$ep[,3]) + 1;
  #   cep.new[,3] = seqlen - rev(ep.new$ep[,2]) + 1;
  #   cep.new[,4] = seqlen - rev(ep.new$ep[,1]) + 1;
  # }      
  
  n = ncol(resmat);
  outcase.exon = c();
  outcase.exon.loc = matrix(0,ncol=n,nrow=nrow(cep.new));
  outcase.exon.sum = outcase.exon.iqr = outcase.exon.loc;
  for (exonnum in 1:nrow(cep.new)) {
    pexon = c(cep.new[exonnum,2]:cep.new[exonnum,3])
    pexon.sum = apply(resmat[pexon,],2,sum)
    pexon.q = quantile(pexon.sum);
    iqr = pexon.q[4]-pexon.q[2];
    iqr.rate = rep(0,n);
    
    lowbound = pexon.q[2]-iqr.cutoff*iqr; upbound = pexon.q[4]+iqr.cutoff*iqr; 
    on.out = which((pexon.sum>upbound) & (apply(GeneData[pexon,],2,FUN=function(x){quantile(x,probs=0.9)})>5)) ;
    off.out = which(pexon.sum<lowbound) ;
    iqr.rate[on.out] = (pexon.sum[on.out] - pexon.q[4])/iqr;
    iqr.rate[off.out] = (pexon.q[2] - pexon.sum[off.out])/iqr;
    
    outcase.tmp = c(on.out,off.out);
    outcase.exon.loc[exonnum,outcase.tmp] =  1;
    outcase.exon.sum[exonnum,] = pexon.sum;
    outcase.exon.iqr[exonnum,] = iqr.rate;
    outcase.exon = unique(c(outcase.exon,outcase.tmp))
    rm(outcase.tmp,iqr.rate);
  }      
  outcase.intron = c();
  outcase.intron.loc = matrix(0,ncol=ncol(resmat),nrow=(nrow(cep.new)-1));
  outcase.intron.sum = outcase.intron.iqr = outcase.intron.loc;
  for (intronnum in 1:(nrow(cep.new)-1)) {
    pintron = c((cep.new[intronnum,3]+5):c(cep.new[(intronnum+1),2]-5)) # from +-2
    pintron.sum = apply(resmat[pintron,],2,sum)
    #         pintron.q = quantile(pintron.sum[which(pintron.sum>0)]);
    pintron.q = quantile(pintron.sum);
    iqr = pintron.q[4]-pintron.q[2];
    iqr.rate = rep(0,n);
    
    lowbound = pintron.q[2]-iqr.cutoff*iqr; upbound = pintron.q[4]+iqr.cutoff*iqr; 
    on.out = which((pintron.sum>upbound) & (apply(GeneData[pintron,],2,FUN=function(x){quantile(x,probs=0.9)})>5)) ;
    off.out = which(pintron.sum<lowbound) ;
    iqr.rate[on.out] = (pintron.sum[on.out] - pintron.q[4])/iqr;
    iqr.rate[off.out] = (pintron.q[2] - pintron.sum[off.out])/iqr;
    
    outcase.tmp = c(on.out,off.out);
    outcase.intron.loc[intronnum,outcase.tmp] =  1;
    outcase.intron.sum[intronnum,] = pintron.sum;
    outcase.intron.iqr[intronnum,] = iqr.rate;
    outcase.intron = unique(c(outcase.intron,outcase.tmp))
    rm(outcase.tmp,iqr.rate);
  }
  outcase = unique(c(outcase.exon,outcase.intron))
  
  return(list(outcase=outcase,cep.new=cep.new,
              outcase.exon=outcase.exon,outcase.intron=outcase.intron,
              outcase.exon.loc=outcase.exon.loc,outcase.intron.loc=outcase.intron.loc,
              outcase.exon.sum=outcase.exon.sum,outcase.intron.sum=outcase.intron.sum,
              outcase.exon.iqr=outcase.exon.iqr,outcase.intron.iqr=outcase.intron.iqr));
}# end function
