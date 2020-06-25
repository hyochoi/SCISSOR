#' Detect local shape changes
#'
#' This function discovers outlying subjects whose RNA-seq have "local" abnormal shapes
#' and provides the most outlying window-direction for each outlier.
#'
#' @param miscGlobalResult Result from \link{miscGlobal}
#' @param pileupData raw coverage matrix, dataI from \link{process_data}
#' @param residualData a residual matrix from PCA, a data matrix subtracted by a low-rank matrix.
#' @param exonset data annotation matrix, exonset from \link{process_data}
#' @param windowSize a window length. Default is 100.
#' @param siglev the significance level for a robust outlier detection. Default is 1e-5.
#' IF cutoff is specified, siglev is not used.
#' @param cutoff the cutoff value for outlying statistics.
#' If NULL, the cutoff value is computed based on the specified siglev.
#' @param reducedReturn
#'
#' @import zoo
#' @export
miscLocal = function(miscGlobalResult,pileupData,Ranges,
                     siglev=1e-4,cutoff=NULL,
                     ADcutoff=3,windowSize=100,
                     reducedReturn=TRUE) {
  n = ncol(miscGlobalResult$residualData); d = nrow(miscGlobalResult$residualData);
  exonset = Ranges$lRanges
  if (nrow(exonset)==1) {
    exon.base=c(exonset[1,2]:exonset[1,3])
  } else {
    exon.base=c()
    for (i in 1:(nrow(exonset)-1)) {
      exon.base=c(exon.base,c(exonset[i,2]:exonset[i,3]))
    }
  }

  out1 = miscGlobalResult$SC;
  residualData2 = miscGlobalResult$residualData[,which(! c(1:n) %in% out1)]
  pileupData2 = pileupData[,which(! c(1:n) %in% out1)]
  MOD=matrix(0,nrow=d,ncol=n)
  NPS=matrix(0,nrow=n,ncol=n)
  OS = rep(0,n);
  if (miscGlobalResult$K==0) {
    out2=out2.sort=c()
    out2stat=rep(0,ncol(residualData2))
    cutoff=0
    residualData_out=NULL
  } else {
    residualData_out = detect_localout(residualData=residualData2,
                                       pileupData=pileupData2,
                                       exonset=exonset,
                                       siglev=siglev,
                                       cutoff=cutoff,
                                       ADcutoff=ADcutoff,
                                       windowSize=windowSize,
                                       reducedReturn=reducedReturn)

    changeID = Relist.hy(1:n,out1);
    out2 = changeID$new[residualData_out$outliers];
    out2.sort = changeID$new[residualData_out$outliers.sort]

    ## Step 2 MOD
    if (length(out1)>=1) {
      MOD[,-out1] = residualData_out$MOD;
      NPS[-out1,-out1] = residualData_out$NPS;
    } else {
      MOD = residualData_out$MOD;
      NPS = NPS
    }
    cutoff = residualData_out$cutoff;
    OS[which(! c(1:n) %in% out1)] = residualData_out$onoff_stat;
  }

  outputObject = list(SC=out2.sort,OS=OS,
                      MOD=MOD,NPS=NPS,
                      cutoff=cutoff,
                      ignoredCases=out1,
                      ks.df=residualData_out$ks.df,ks.pval=residualData_out$ks.pval,
                      onoffResults=residualData_out)
  class(outputObject) = append(class(outputObject),"LSCOutput")
  return(outputObject)
}

#' Detect local shape outliers with a residual matrix
#'
#' This function identifies local shape variants by mining locally on/off abnormalities.
#'
#'
#' @export
detect_localout = function(residualData,pileupData,exonset,
                           siglev=NULL,cutoff=NULL,ADcutoff=3,
                           windowSize=100,
                           reducedReturn=FALSE) {

  n2=ncol(residualData)
  onoff_res=as.list(NULL)
  readconstr=10 # the minimum reads count required to be considered for on/off local shape changes.
  onoff_res[[1]]=get_offstat(residualData=residualData,pileupData=pileupData,exonset=exonset,ADcutoff=ADcutoff,
                             windowSize=windowSize,readconstr=readconstr)
  onoff_res[[2]]=get_exon_offstat(residualData=residualData,pileupData=pileupData,exonset=exonset,ADcutoff=ADcutoff)
  onoff_res[[3]]=get_exon_onstat(residualData=residualData,pileupData=pileupData,
                                 exonset=exonset,readconstr=readconstr)
  onoff_res[[4]]=get_intron_onstat(residualData=residualData,pileupData=pileupData,
                                   exonset=exonset,readconstr=readconstr)

  onoff_statmat=matrix(0,nrow=length(onoff_res),ncol=n2)
  for (i in 1:length(onoff_res)) {
    onoff_statmat[i,]=onoff_res[[i]]$stat
  }
  onoff_stat=apply(onoff_statmat,2,max)
  onoff_stat_which=apply(onoff_statmat,2,which.max)

  MOD = matrix(0,nrow=nrow(residualData),ncol=n2)
  NPS = matrix(0,nrow=n2,ncol=n2)
  for (j in 1:n2) {
    MOD[,j]=onoff_res[[onoff_stat_which[j]]]$tMOD[,j]
    NPS[,j]=onoff_res[[onoff_stat_which[j]]]$tNPS[,j]
  }

  ks.pval=ks.stat=rep(0,20)
  for (df in 1:20) {
    # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
    temp.out=which(onoff_stat>sqrt(qchisq(siglev,df=df,lower.tail=F)))
    ks.output=ks.test(onoff_stat[which(! c(1:n2) %in% temp.out)]^2,pchisq,df)
    ks.pval[df]=ks.output$p.value
    ks.stat[df]=ks.output$statistic
    # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
    #           round(ks.stat[df],digits=3)),"\n")
  }

  df=which.min(ks.stat)
  if (is.null(cutoff)) {
    cutoff=sqrt(qchisq(siglev,df=df,lower.tail=F))
  }
  onoffout=c()
  for (i in 1:length(onoff_res)) {
    onoffout=unique(c(onoffout,which(onoff_res[[i]]$stat>cutoff)))
  }
  if (length(onoffout)>0) {
    onoffout.sort=onoffout[order(onoff_stat[onoffout],decreasing=T)]
  } else {
    onoffout=NULL; onoffout.sort=NULL;
  }
  if (reducedReturn) {
    onoff_res=NULL
  }
  return(list(outliers=onoffout,
              outliers.sort=onoffout.sort,
              onoff_stat=onoff_stat,
              MOD=MOD,NPS=NPS,
              cutoff=cutoff,ks.df=df,ks.pval=ks.pval[df],
              onoff_res=onoff_res))
}

#' Quantifying abnormal deletions in narrow regions
#'
#' This function quantifies abnormality associated with
#' @export
get_offstat= function(residualData,pileupData,exonset,ADcutoff=3,
                      windowSize=100,readconstr=10) {
  ## Step 2 : detect "off" outliers
  require(zoo)
  # Find eligible regions
  n2 = ncol(residualData)
  medvec = apply(pileupData,1,median)

  if (nrow(exonset)==1) {
    exon.base=c(exonset[1,2]:exonset[1,3])
  } else {
    exon.base=c()
    for (i in 1:(nrow(exonset)-1)) {
      exon.base=c(exon.base,c(exonset[i,2]:exonset[i,3]))
    }
  }
  std.overall=t(apply(residualData[exon.base,],1,pd.rate.hy))
  mad.overall=apply(std.overall,2,mad)
  mad.overall[which(mad.overall<1)]=1
  residualData = sweep(residualData,2,mad.overall,"/");

  medvec.case = apply(pileupData[exon.base,],2,median)
  medvec.case01 = rep(0,n2); medvec.case01[which(medvec.case>10)]=1
  window_min = rollapply(medvec,width=windowSize,FUN=min) # min of med
  # window_mean = rollapply(medvec,width=windowSize,FUN=mean)

  window_idx0 = which(window_min>=readconstr)
  window_idx1 = window_idx0[which(window_idx0 %in% exon.base)]
  window_idx = window_idx1[which(window_min[window_idx1]>=quantile(window_min[window_idx1],probs=0.1))] # do we really need this?

  # Find exons whose lengths are less than the pre-determined size.
  smallexon_0 = which((exonset[1:(nrow(exonset)-1),3]-exonset[1:(nrow(exonset)-1),2]+1)<windowSize)
  smallexon_min = c()
  for (j in smallexon_0) {
    smallexon_min = min(medvec[c(exonset[j,2]:exonset[j,3])])
  }
  smallexon = smallexon_0[which(smallexon_min>=readconstr)]

  # Determine window size
  winsize_vec = rep(0,length(window_min))
  winsize_vec[window_idx] = windowSize;
  winsize_vec[exonset[smallexon,2]] = exonset[smallexon,3]-exonset[smallexon,2]+1

  # Final start_positions for calculating statistics
  window_site = which(winsize_vec>0)
  window_size = winsize_vec[window_site]
  pdrate = matrix(0,nrow=length(window_site),ncol=ncol(residualData))
  tMOD=matrix(0,nrow=nrow(pileupData),ncol=n2) # temp MOD
  tNPS=matrix(0,nrow=n2,ncol=n2) # temp NPS
  if ((length(window_site)>0) & (length(which(medvec.case01==1))>30)) {
    for (i in which((window_site>100) & (window_site<(nrow(residualData)-100)))) {
      carea = c((window_site[i]):(window_site[i]-1+window_size[i]))
      # pdrate[i,] = pd.rate.hy(-apply(residualData[carea,],2,sum),qrsc=T)
      pdrate[i,which(medvec.case01==1)] = pd.rate.hy(-apply(residualData[carea,which(medvec.case01==1)],2,sum),qrsc=F)
    }
    outdir = which(apply(pdrate[,which(medvec.case01==1)],1,ADstatWins.hy)>=ADcutoff)
    pdrate[outdir,]=0

    # pdrate2 = sweep(pdrate,2,medvec.case01,"*")
    pdrate3 = pdrate
    off_stat = apply(pdrate3,2,max)
    where_on = apply(pdrate3,2,which.max)

    for (j in 1:n2) {
      start_loc = window_site[where_on[j]]
      end_loc = start_loc + window_size[where_on[j]] - 1
      carea = c(start_loc:end_loc)
      tMOD[carea,j] = -1/sqrt(length(carea))
      tNPS[,j] = pdrate3[where_on[j],]
    }
  } else {
    off_stat = rep(0,ncol(residualData))
  }

  return(list(stat=off_stat,tMOD=tMOD,tNPS=tNPS))
}

#'
#' @export
get_exon_offstat= function(residualData,pileupData,exonset,ADcutoff=3) {

  n2 = ncol(residualData)
  medvec = apply(pileupData,1,median)
  nexon=nrow(exonset)

  exon.comb=diag(1,nexon)
  if (nexon>2) {
    for (j in 2:(nexon-1)) {
      temp=rep(0,nexon)
      temp[c(1:j)]=1
      exon.comb=cbind(exon.comb,temp,rev(temp))
    }
  }
  colnames(exon.comb)=NULL

  pdrate = matrix(0,nrow=ncol(exon.comb),ncol=n2)
  for (i in 1:nrow(pdrate)) {
    carea=c()
    for (e in which(exon.comb[,i]==1)) {
      carea=c(carea,c(exonset[e,2]:exonset[e,3]))
    }
    if (length(carea)>20) {
      pdrate[i,] = pd.rate.hy(-apply(residualData[carea,],2,sum),qrsc=F);
    }
  }

  outdir=which(apply(pdrate,1,ADstatWins.hy)>=ADcutoff)
  pdrate[outdir,]=0
  exonoff_stat = apply(pdrate,2,max)
  where_on = apply(pdrate,2,which.max)

  tMOD=matrix(0,nrow=nrow(pileupData),ncol=n2) # temp MOD
  tNPS=matrix(0,nrow=n2,ncol=n2) # temp NPS
  for (j in 1:n2) {
    carea=c()
    for (e in which(exon.comb[,where_on[j]]==1)) {
      carea=c(carea,c(exonset[e,2]:exonset[e,3]))
    }
    if (length(carea)>20) {
      tMOD[carea,j] = -1/sqrt(length(carea))
      tNPS[,j] = pdrate[where_on[j],]
    }
  }

  return(list(stat=exonoff_stat,tMOD=tMOD,tNPS=tNPS))
}

#'
#' @export
get_exon_onstat= function(residualData,pileupData,exonset,readconstr=10) {

  n2 = ncol(residualData)
  medvec = apply(pileupData,1,median)
  nexon=nrow(exonset)

  exon.comb=diag(1,nexon)
  if (nexon>2) {
    for (j in 2:(nexon-1)) {
      temp=rep(0,nexon)
      temp[c(1:j)]=1
      exon.comb=cbind(exon.comb,temp,rev(temp))
    }
  }
  colnames(exon.comb)=NULL

  medianmat = matrix(0,nrow=ncol(exon.comb),ncol=n2)
  medianmat1 = matrix(0,nrow=ncol(exon.comb),ncol=n2)
  medianmat2 = matrix(0,nrow=ncol(exon.comb),ncol=n2)
  pdrate1 = matrix(0,nrow=ncol(exon.comb),ncol=n2)
  pdrate2 = matrix(0,nrow=ncol(exon.comb),ncol=n2)
  for (i in 1:ncol(exon.comb)) {

    carea1=c()
    for (e in which(exon.comb[,i]==1)) {
      carea1=c(carea1,c(exonset[e,2]:exonset[e,3]))
    }
    if (length(carea1)>20) {
      pdrate1[i,] = pd.rate.hy(apply(residualData[carea1,],2,sum),qrsc=T);
      medianmat[i,] = apply(pileupData[carea1,],2,FUN=function(x){(quantile(x,probs=0.75)>=readconstr)})
      medianmat1[i,] = apply(pileupData[carea1,],2,FUN=function(x){fivenum(x)[3]})
    }

    # carea2=c()
    # for (e in which(exon.comb[,i]==0)) {
    #   carea2=c(carea2,c(exonset[e,2]:exonset[e,3]))
    # }
    # pdrate2[i,] = pd.rate.hy(apply(residualData[carea2,],2,sum),qrsc=T);
    # medianmat2[i,] = apply(pileupData[carea2,],2,FUN=function(x){fivenum(x)[3]})
  }

  exonon_stat = apply(pdrate1*medianmat,2,max)
  where_on = apply(pdrate1*medianmat,2,which.max)

  tMOD=matrix(0,nrow=nrow(pileupData),ncol=n2)
  tNPS=matrix(0,nrow=n2,ncol=n2) # temp NPS
  for (j in 1:n2) {
    carea=c()
    for (e in which(exon.comb[,where_on[j]]==1)) {
      carea=c(carea,c(exonset[e,2]:exonset[e,3]))
    }
    if (length(carea)>20) {
      tMOD[carea,j] = 1/sqrt(length(carea))
      tNPS[,j] = pdrate1[where_on[j],]
    }
  }

  return(list(stat=exonon_stat,tMOD=tMOD,tNPS=tNPS))
}

#'
#' @export
get_intron_onstat= function(residualData,pileupData,exonset,readconstr=10) {

  intron_onstat_inner2=function(x){
    lcarea=length(carea)
    y=max(quantile(x[carea[6:ceiling(lcarea/3)]],probs=0.75),
          quantile(x[carea[ceiling(lcarea/3):(2*ceiling(lcarea/3))]],probs=0.75),
          quantile(x[carea[(2*ceiling(lcarea/3)):(lcarea-6)]],probs=0.75))
    max_exonread = max(quantile(x[carea1],probs=0.5),quantile(x[carea2],probs=0.5));
    if (max_exonread>readconstr) {
      readconstr2=max(5,0.1*max_exonread)
    } else {
      readconstr2=readconstr
    }
    if (max_exonread>100) {
      ((y>=readconstr) & (y>=readconstr2))
    } else {
      ((y>=readconstr) | (y>=readconstr2))
    }
  }

  n2 = ncol(residualData)
  nexon = nrow(exonset)

  medianmat = matrix(0,nrow=(nexon-1),ncol=n2)
  pdrate = matrix(0,nrow=(nexon-1),ncol=n2)
  if (nexon > 1) {
    for (i in 1:(nexon-1)) {

      carea = c((exonset[i,3]+1):(exonset[(i+1),2]-1))
      if (length(carea)>20) {
        carea1 = c(exonset[i,2]:exonset[i,3])
        carea2 = c(exonset[(i+1),2]:exonset[(i+1),3])
        medianmat[i,] = apply(pileupData,2,intron_onstat_inner2)

        pdrate[i,] = pd.rate.hy(apply(residualData[carea,],2,sum),qrsc=T);
      }
    }
  }

  intron_onstat = apply(pdrate*medianmat,2,max)
  where_on = apply(pdrate*medianmat,2,which.max)

  tMOD=matrix(0,nrow=nrow(pileupData),ncol=n2)
  tNPS=matrix(0,nrow=n2,ncol=n2) # temp NPS
  if (nexon > 1) {
    for (j in 1:n2) {
      carea=c((exonset[where_on[j],3]+1):(exonset[(where_on[j]+1),2]-1))
      tMOD[carea,j] = 1/sqrt(length(carea))
      tNPS[,j] = pdrate[where_on[j],]
    }
  }

  return(list(stat=intron_onstat,tMOD=tMOD,tNPS=tNPS))
}
