
miscGlobal_test5 = function(inputData,pileupData,Ranges,JSR.table,AToutlier=T,
                            siglev=1e-04,ADcutoff=2,print.proc=F) {

  n = dim(inputData)[2]
  plat.table = build_blockTable(Ranges=Ranges)
  plat.baseMat = build_baseMat(plat.table)
  size.adjMat = sqrt(t(plat.baseMat)%*%plat.baseMat)

  # Estimate exons
  exonlabels = label_exons(pileupData=pileupData,Ranges=Ranges)
  estexons = names(which(exonlabels==1))

  knownDir = build_knownDir(plat.table=plat.table,Ranges=Ranges,JSR.table=JSR.table,exons=estexons)
  adj.knownDir = size.adjMat%*%knownDir # adjust size of plats
  plat.names = rownames(knownDir)
  plat.names2 = sapply(plat.names,FUN=function(t){unlist(strsplit(t,"[.]"))[1]})
  platBasisDir = apply(plat.baseMat,2,get_unitdir)
  KnownBasisDir = apply(plat.baseMat%*%knownDir,2,get_unitdir)
  normProjData = t(platBasisDir) %*% inputData # low-dimensional data object to be used for outlier detection
  # normProjData = t(apply(normProjData,1,make_normal))

  ## 1. Detecting outliers from known directions
  ## 1-1. Get projection outlyingness
  # Get primary projection outlyingness
  knownPOmat = get_POgivenB2(X=normProjData,B=adj.knownDir,qrsc=TRUE,ADcutoff=ADcutoff)
  knownPO = diag(knownPOmat[apply(abs(knownPOmat),2,which.max),])
  # knownPODir = knownDir[,apply(abs(knownPOmat),2,which.max)]

  ## 1-2. Detecting non-outliers
  intron.nonout = get_intron_nonout(pileupData=pileupData,Range=Ranges)
  exon.nonout = get_exon_nonout(normData=inputData,POmat=knownPOmat,Ranges=Ranges,plat.table=plat.table)
  ## 1-3. Non-unimodality / Extreme skewness
  # crypPOmat = sweep(x=knownPOmat,1,STATS=apply(crypPOmat,1,filter_multimodal),"*")

  ## 2. Get sum of squares of the residuals
  knownROmat = get_resdZsum(X=normProjData,B=adj.knownDir,L1=FALSE,qrsc=TRUE) # Residual (Z-score) outlyingness
  knownROmat0 = knownROmat
  for (ii in which(grepl("I",rownames(knownROmat)))) {
    tmp.name = unlist(strsplit(rownames(knownROmat)[ii],"[.]"))[1]
    knownROmat0[ii,which(intron.nonout[which(rownames(intron.nonout)==tmp.name),]==0)] = 1000
  }
  for (ii in which(grepl("E",rownames(knownROmat)))) {
    tmp.name = unlist(strsplit(rownames(knownROmat)[ii],"[.]"))[1]
    knownROmat0[ii,which(exon.nonout[which(rownames(exon.nonout)==tmp.name),]==0)] = 1000
  }
  for (ii in which(grepl("SV",rownames(knownROmat)))) {
    tmp.region = rownames(knownDir)[which(! knownDir[,ii]==0)]
    tmp.case = which(apply(matrix(exon.nonout[which(rownames(exon.nonout) %in% tmp.region),],nrow=length(tmp.region)),2,function(x) length(which(x==0)))>0)
    knownROmat0[ii,tmp.case] = 1000
  }
  for (ii in which(grepl("J",rownames(knownROmat)))) {
    tmp.region = rownames(knownDir)[which(! knownDir[,ii]==0)]
    tmp.case = which(apply(matrix(exon.nonout[which(rownames(exon.nonout) %in% tmp.region),],nrow=length(tmp.region)),2,function(x) length(which(x==0)))>0)
    knownROmat0[ii,tmp.case] = 1000
  }

  knownRO0 = diag(knownROmat[apply(knownROmat0,2,which.min),])
  knownRODir0 = knownDir[,apply(knownROmat0,2,which.min)]
  # knownPOminRO0 = diag(knownPOmat[apply(knownROmat,2,which.min),])
  # knownPOminRO0 = diag(get_POgivenB2(X=normProjData,B=size.adjMat%*%knownRODir0,qrsc=TRUE,ADcutoff=ADcutoff))
  knownRODir = knownRODir0

  X.stats = apply(normProjData,1,FUN=function(t){unlist(pd.stats(t,qrsc=TRUE))})
  X.stats[3,which(grepl(pattern="I",rownames(normProjData))==T)] = X.stats[2,which(grepl(pattern="I",rownames(normProjData))==T)]
  ADstat = apply(normProjData,1,ADstatWins.hy)
  ADw = sapply(ADstat,FUN=function(t){max(1-(t/100),0.1)})
  st = Sys.time()
  for (case in 1:n) {
    tmpdir0 = knownRODir0[,case]
    tmpRO0 = knownRO0[case]
    # tmpPO0 = knownPOminRO0[case]
    init.min = min(which(tmpdir0==1))
    init.max = max(which(tmpdir0==1))

    if (init.min>1) {
      for (i in (init.min-1):1) {
        if (plat.names2[i] %in% rownames(intron.nonout)[which(intron.nonout[,case]==0)]) {
          next
        }
        if (plat.names2[i] %in% rownames(exon.nonout)[which(exon.nonout[,case]==0)]) {
          next
        }
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        resdX = normProjData - direction1%*%t(direction1)%*%normProjData
        tmpRO1sq = 0
        for (j in 1:nrow(resdX)) {
          tmpRO1sq = tmpRO1sq + ADw[j]*(sapply(resdX[j,case],FUN=function(t) {compute.pd(x=t,stats=X.stats[,j])})^2)
        }
        tmpRO1 = sqrt(tmpRO1sq)
        # tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,cases=case,qrsc=TRUE)
        # tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.0001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          # tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }

    if (init.max<length(tmpdir0)) {
      for (i in (init.max+1):length(tmpdir0)) {
        if (plat.names2[i] %in% rownames(intron.nonout)[which(intron.nonout[,case]==0)]) {
          next
        }
        if (plat.names2[i] %in% rownames(exon.nonout)[which(exon.nonout[,case]==0)]) {
          next
        }
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        resdX = normProjData - direction1%*%t(direction1)%*%normProjData
        tmpRO1sq = 0
        for (j in 1:nrow(resdX)) {
          tmpRO1sq = tmpRO1sq + ADw[j]*(sapply(resdX[j,case],FUN=function(t) {compute.pd(x=t,stats=X.stats[,j])})^2)
        }
        tmpRO1 = sqrt(tmpRO1sq)
        # tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,cases=case,qrsc=TRUE)
        # tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.0001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          # tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }
    if (sum((knownRODir0[,case]-tmpdir0)^2)>0) {
      knownRODir[,case] = tmpdir0
      # colnames(knownRODir)[case] = strsplit(rownames(knownRODir)[min(which(knownRODir[,case]==1))],"[.]")[[1]][1]
    }
    rm(tmpdir0)
    if (print.proc) {
      cat(case,"\n")
    }
  }
  et = Sys.time()
  et - st
  knownPOminRO = diag(get_POgivenB2(X=normProjData,B=size.adjMat%*%knownRODir,qrsc=TRUE,ADcutoff=ADcutoff))
  knownRO = diag(get_resdZsum(X=normProjData,B=size.adjMat%*%knownRODir,L1=FALSE)) # Residual (Z-score) outlyingness
  # knownPOminRO = diag(knownPOmat[apply(knownROmat,2,which.min),])

  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(knownPOminRO^2 < 1e-10)
  for (df in 1:20) {
    # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
    temp.out=which(abs(knownPOminRO)>sqrt(qchisq(0.001,df=df,lower.tail=F)))
    ks.output=ks.test(unique((knownPOminRO^2)[which(! c(1:length(knownPOminRO)) %in% c(temp.out,temp.zero))]),pchisq,df)
    ks.pval[df]=ks.output$p.value
    ks.stat[df]=ks.output$statistic
    # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
    #           round(ks.stat[df],digits=3)),"\n")
  }
  df=which.min(ks.stat)
  cutoff = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),5)
  # cutoff = sqrt(qchisq(p=(1-siglev),df=df))

  # Save results
  OStmp = knownPOminRO
  MODtmp = knownRODir
  OUTtmp = which(abs(OStmp)>cutoff)

  ## 2. Detecting outliers from ATS/ATT directions
  if (AToutlier) {
    atDir = build_atDir(plat.table=plat.table,
                        inputData=normData,Ranges=Ranges,JSR.table=JSR.table,
                        cutoff=cutoff)
    if (ncol(atDir)>0) {
      adj.atDir = size.adjMat%*%atDir
      # atPOmat = get_POgivenB(X=normProjData,B=adj.atDir,qrsc=TRUE)
      # atPO = diag(atPOmat[apply(abs(atPOmat),2,which.max),])
      # atPODir = atDir[,apply(atPOmat,2,which.max)]
      atROmat = get_resdZsum(X=normProjData,B=adj.atDir,L1=FALSE,qrsc=TRUE)
      atROmat0 = atROmat
      for (ii in 1:nrow(atROmat)) {
        tmp.name = gsub("ATS","I",unlist(strsplit(rownames(atROmat)[ii],"[.]"))[1])
        tmp.name = gsub("ATT","I",tmp.name)
        atROmat0[ii,which(intron.nonout[which(rownames(intron.nonout)==tmp.name),]==0)] = 1000
      }
      atRO = diag(atROmat[apply(atROmat0,2,which.min),])
      atRODir = atDir[,apply(atROmat0,2,which.min)]
      atPOminRO = diag(get_POgivenB2(X=normProjData,B=size.adjMat%*%atRODir,qrsc=TRUE,ADcutoff=ADcutoff))
      # atPOminRO = diag(atPOmat[apply(atROmat0,2,which.min),])
      # atSC = which((abs(atPOminRO)>cutoff) & ((abs(atPOminRO)>abs(OStmp)) & (atRO<knownRO)))
      atSC = which((atPOminRO>cutoff) & (atRO<knownRO))
      atSC = atSC[which((knownRO^2 - atRO^2)[atSC]>qchisq(p=0.0001,df=1,lower.tail=F))]

      ## Update the OS and MOD
      OStmp[atSC] = atPOminRO[atSC]
      MODtmp[,atSC] = atRODir[,atSC]
      colnames(MODtmp)[atSC] = colnames(atRODir)[atSC]
    }
  }

  ## Get final results table
  # Find region tags
  cite_first = function(t) {
    return(strsplit(t,"[.]")[[1]][1])
  }
  find_region.tag = function(MOD) {
    region.names = rownames(MOD)
    inner_fn = function(x) {
      paste(unique(sapply(region.names[which(! x==0)],cite_first)),collapse=",")
    }
    # region.tag1 = sapply(colnames(MOD),cite_first)
    # region.tag2 = apply(MOD,2,inner_fn)
    # region.tag0 = paste(region.tag1,region.tag2,sep=",")
    # region.tag = sapply(region.tag0,FUN=function(t){paste(unique(unlist(strsplit(t,","))),collapse=",")})
    region.tag = apply(MOD,2,inner_fn)
    return(region.tag)
  }
  region.tag = find_region.tag(MOD=MODtmp)
  # paste(unique(sapply(region.names[which(! MODtmp[,4]==0)],cite_first)),collapse=",")

  # Results
  outliers = which(abs(OStmp)>cutoff)
  output.tag = region.tag[outliers]
  basis.tag = sapply(colnames(MODtmp)[outliers],cite_first)
  # output.tag = colnames(MODtmp[,outliers])
  # if (length(which(grepl("[.]",output.tag)))>0) {
  #   output.tag[which(grepl("[.]",output.tag))] =
  #     sapply(output.tag[which(grepl("[.]",output.tag))],FUN=function(t){strsplit(t,"[.]")[[1]][1]})
  # }
  globalResult.table = data.frame(Outlier=outliers,
                                  Statistic=round(abs(OStmp[outliers]),digits=3),
                                  Region=locate_region(tag=colnames(MODtmp)[outliers],Ranges=Ranges,JSR.table=JSR.table),
                                  Region.tag=output.tag,
                                  Basis.tag=basis.tag)

  ## Output
  outputObject = list(table=globalResult.table,
                      Outlier=globalResult.table$Outlier,
                      OS=abs(OStmp),
                      signOS=OStmp,
                      platBasis=plat.baseMat,
                      platMOD=MODtmp,
                      cutoff=cutoff)
  class(outputObject) = append(class(outputObject),"GSCOutput")
  return(outputObject)
}

miscGlobal_test4 = function(inputData,pileupData,Ranges,JSR.table,
                            siglev=1e-04,ADcutoff=4) {

  n = dim(inputData)[2]
  plat.table = build_blockTable(Ranges=Ranges)
  plat.baseMat = build_baseMat(plat.table)
  size.adjMat = sqrt(t(plat.baseMat)%*%plat.baseMat)

  knownDir = build_knownDir(plat.table=plat.table,Ranges=Ranges,JSR.table=JSR.table)
  adj.knownDir = size.adjMat%*%knownDir # adjust size of plats
  plat.names = rownames(knownDir)
  plat.names2 = sapply(plat.names,FUN=function(t){unlist(strsplit(t,"[.]"))[1]})
  platBasisDir = apply(plat.baseMat,2,get_unitdir)
  KnownBasisDir = apply(plat.baseMat%*%knownDir,2,get_unitdir)
  normProjData = t(platBasisDir) %*% inputData # low-dimensional data object to be used for outlier detection
  # normProjData = t(apply(normProjData,1,make_normal))

  ## 0. Detecting non-outliers from intronic regions
  intron.nonout = get_intron_nonout(pileupData=pileupData,Range=Ranges)
  ## 1. Detecting outliers from known directions

  ## 1-1. Get projection outlyingness
  # Get primary projection outlyingness
  knownPOmat = get_POgivenB(X=normProjData,B=adj.knownDir,qrsc=TRUE)
  knownPO = diag(knownPOmat[apply(abs(knownPOmat),2,which.max),])
  # knownPODir = knownDir[,apply(abs(knownPOmat),2,which.max)]

  ## 1-2. Non-unimodality / Extreme skewness
  # crypPOmat = sweep(x=knownPOmat,1,STATS=apply(crypPOmat,1,filter_multimodal),"*")

  ## 2. Get sum of squares of the residuals
  knownROmat = get_resdZsum(X=normProjData,B=adj.knownDir,L1=FALSE,qrsc=TRUE) # Residual (Z-score) outlyingness
  knownROmat0 = knownROmat
  for (ii in which(grepl("I",rownames(knownROmat)))) {
    tmp.name = unlist(strsplit(rownames(knownROmat)[ii],"[.]"))[1]
    knownROmat0[ii,which(intron.nonout[which(rownames(intron.nonout)==tmp.name),]==0)] = 1000
  }
  knownRO0 = diag(knownROmat[apply(knownROmat0,2,which.min),])
  knownRODir0 = knownDir[,apply(knownROmat0,2,which.min)]
  # knownPOminRO0 = diag(knownPOmat[apply(knownROmat,2,which.min),])
  # knownPOminRO0 = diag(get_POgivenB2(X=normProjData,B=size.adjMat%*%knownRODir0,qrsc=TRUE,ADcutoff=ADcutoff))
  knownRODir = knownRODir0

  for (case in 1:n) {
    tmpdir0 = knownRODir0[,case]
    tmpRO0 = knownRO0[case]
    # tmpPO0 = knownPOminRO0[case]
    init.min = min(which(tmpdir0==1))
    init.max = max(which(tmpdir0==1))

    if (init.min>1) {
      for (i in (init.min-1):1) {
        if (plat.names2[i] %in% rownames(intron.nonout)[which(intron.nonout[,case]==0)]) {
          next
        }
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,qrsc=TRUE)[case]
        # tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.0001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          # tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }

    if (init.max<length(tmpdir0)) {
      for (i in (init.max+1):length(tmpdir0)) {
        if (plat.names2[i] %in% rownames(intron.nonout)[which(intron.nonout[,case]==0)]) {
          next
        }
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,qrsc=TRUE)[case]
        # tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.0001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          # tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }
    if (sum((knownRODir0[,case]-tmpdir0)^2)>0) {
      knownRODir[,case] = tmpdir0
      # colnames(knownRODir)[case] = strsplit(rownames(knownRODir)[min(which(knownRODir[,case]==1))],"[.]")[[1]][1]
    }
    rm(tmpdir0)
    cat(case,"\n")
  }
  knownPOminRO = diag(get_POgivenB2(X=normProjData,B=size.adjMat%*%knownRODir,qrsc=TRUE,ADcutoff=ADcutoff))
  knownRO = diag(get_resdZsum(X=normProjData,B=size.adjMat%*%knownRODir,L1=FALSE)) # Residual (Z-score) outlyingness
  # knownPOminRO = diag(knownPOmat[apply(knownROmat,2,which.min),])

  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(knownPOminRO^2 < 1e-10)
  for (df in 1:20) {
    # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
    temp.out=which(abs(knownPOminRO)>sqrt(qchisq(0.001,df=df,lower.tail=F)))
    ks.output=ks.test(unique((knownPOminRO^2)[which(! c(1:length(knownPOminRO)) %in% c(temp.out,temp.zero))]),pchisq,df)
    ks.pval[df]=ks.output$p.value
    ks.stat[df]=ks.output$statistic
    # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
    #           round(ks.stat[df],digits=3)),"\n")
  }
  df=which.min(ks.stat)
  cutoff = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),4)
  # cutoff = sqrt(qchisq(p=(1-siglev),df=df))

  # Save results
  OStmp = knownPOminRO
  MODtmp = knownRODir
  OUTtmp = which(abs(OStmp)>cutoff)

  ## 2. Detecting outliers from ATS/ATT directions
  atDir = build_atDir(plat.table=plat.table,
                      inputData=normData,Ranges=Ranges,JSR.table=JSR.table,
                      cutoff=cutoff)
  if (ncol(atDir)>0) {
    adj.atDir = size.adjMat%*%atDir
    # atPOmat = get_POgivenB(X=normProjData,B=adj.atDir,qrsc=TRUE)
    # atPO = diag(atPOmat[apply(abs(atPOmat),2,which.max),])
    # atPODir = atDir[,apply(atPOmat,2,which.max)]
    atROmat = get_resdZsum(X=normProjData,B=adj.atDir,L1=FALSE,qrsc=TRUE)
    atROmat0 = atROmat
    for (ii in 1:nrow(atROmat)) {
      tmp.name = gsub("ATS","I",unlist(strsplit(rownames(atROmat)[ii],"[.]"))[1])
      tmp.name = gsub("ATT","I",tmp.name)
      atROmat0[ii,which(intron.nonout[which(rownames(intron.nonout)==tmp.name),]==0)] = 1000
    }
    atRO = diag(atROmat[apply(atROmat0,2,which.min),])
    atRODir = atDir[,apply(atROmat0,2,which.min)]
    atPOminRO = diag(get_POgivenB2(X=normProjData,B=size.adjMat%*%atRODir,qrsc=TRUE,ADcutoff=ADcutoff))
    # atPOminRO = diag(atPOmat[apply(atROmat0,2,which.min),])
    # atSC = which((abs(atPOminRO)>cutoff) & ((abs(atPOminRO)>abs(OStmp)) & (atRO<knownRO)))
    atSC = which((atPOminRO>cutoff) & (atRO<knownRO))
    atSC = atSC[which((knownRO^2 - atRO^2)[atSC]>qchisq(p=0.0001,df=1,lower.tail=F))]

    ## Update the OS and MOD
    OStmp[atSC] = atPOminRO[atSC]
    MODtmp[,atSC] = atRODir[,atSC]
    colnames(MODtmp)[atSC] = colnames(atRODir)[atSC]
  }


  ## Get final results table
  # Find region tags
  cite_first = function(t) {
    return(strsplit(t,"[.]")[[1]][1])
  }
  find_region.tag = function(MOD) {
    region.names = rownames(MOD)
    inner_fn = function(x) {
      paste(unique(sapply(region.names[which(! x==0)],cite_first)),collapse=",")
    }
    # region.tag1 = sapply(colnames(MOD),cite_first)
    # region.tag2 = apply(MOD,2,inner_fn)
    # region.tag0 = paste(region.tag1,region.tag2,sep=",")
    # region.tag = sapply(region.tag0,FUN=function(t){paste(unique(unlist(strsplit(t,","))),collapse=",")})
    region.tag = apply(MOD,2,inner_fn)
    return(region.tag)
  }
  region.tag = find_region.tag(MOD=MODtmp)
  # paste(unique(sapply(region.names[which(! MODtmp[,4]==0)],cite_first)),collapse=",")

  # Results
  outliers = which(abs(OStmp)>cutoff)
  output.tag = region.tag[outliers]
  basis.tag = sapply(colnames(MODtmp)[outliers],cite_first)
  # output.tag = colnames(MODtmp[,outliers])
  # if (length(which(grepl("[.]",output.tag)))>0) {
  #   output.tag[which(grepl("[.]",output.tag))] =
  #     sapply(output.tag[which(grepl("[.]",output.tag))],FUN=function(t){strsplit(t,"[.]")[[1]][1]})
  # }
  globalResult.table = data.frame(Outlier=outliers,
                                  Statistic=round(abs(OStmp[outliers]),digits=3),
                                  Region=locate_region(tag=colnames(MODtmp[,outliers]),Ranges=Ranges,JSR.table=JSR.table),
                                  Region.tag=output.tag,
                                  Basis.tag=basis.tag)

  ## Output
  outputObject = list(table=globalResult.table,
                      Outlier=globalResult.table$Outlier,
                      OS=abs(OStmp),
                      signOS=OStmp,
                      platBasis=plat.baseMat,
                      platMOD=MODtmp,
                      cutoff=cutoff)
  class(outputObject) = append(class(outputObject),"GSCOutput")
  return(outputObject)
}

miscGlobal_test3 = function(inputData,pileupData,Ranges,JSR.table,
                            siglev=1e-04,ADcutoff=3) {

  n = dim(inputData)[2]
  plat.table = build_blockTable(Ranges=Ranges)
  plat.baseMat = build_baseMat(plat.table)
  size.adjMat = sqrt(t(plat.baseMat)%*%plat.baseMat)

  knownDir = build_knownDir(plat.table=plat.table,Ranges=Ranges,JSR.table=JSR.table)
  adj.knownDir = size.adjMat%*%knownDir # adjust size of plats
  plat.names = rownames(knownDir)
  plat.names2 = sapply(plat.names,FUN=function(t){unlist(strsplit(t,"[.]"))[1]})
  platBasisDir = apply(plat.baseMat,2,get_unitdir)
  KnownBasisDir = apply(plat.baseMat%*%knownDir,2,get_unitdir)
  normProjData = t(platBasisDir) %*% normData # low-dimensional data object to be used for outlier detection
  # normProjData = t(apply(normProjData,1,make_normal))

  ## 0. Detecting non-outliers from intronic regions
  intron.nonout = get_intron_nonout(pileupData=rawData,Range=Ranges)
  ## 1. Detecting outliers from known directions

  ## 1-1. Get projection outlyingness
  # Get primary projection outlyingness
  knownPOmat = get_POgivenB(X=normProjData,B=adj.knownDir,qrsc=TRUE)
  knownPO = diag(knownPOmat[apply(abs(knownPOmat),2,which.max),])
  # knownPODir = knownDir[,apply(abs(knownPOmat),2,which.max)]

  ## 1-2. Non-unimodality / Extreme skewness
  # crypPOmat = sweep(x=knownPOmat,1,STATS=apply(crypPOmat,1,filter_multimodal),"*")

  ## 2. Get sum of squares of the residuals
  knownROmat = get_resdZsum(X=normProjData,B=adj.knownDir,L1=FALSE) # Residual (Z-score) outlyingness
  knownROmat0 = knownROmat
  for (ii in which(grepl("I",rownames(knownROmat)))) {
    tmp.name = unlist(strsplit(rownames(knownROmat)[ii],"[.]"))[1]
    knownROmat0[ii,which(intron.nonout[which(rownames(intron.nonout)==tmp.name),]==0)] = 1000
  }
  knownRO0 = diag(knownROmat[apply(knownROmat0,2,which.min),])
  knownRODir0 = knownDir[,apply(knownROmat0,2,which.min)]
  # knownPOminRO0 = diag(knownPOmat[apply(knownROmat,2,which.min),])
  knownPOminRO0 = diag(get_POgivenB(X=normProjData,B=size.adjMat%*%knownRODir0,qrsc=TRUE))
  knownRODir = knownRODir0

  for (case in 1:n) {
    tmpdir0 = knownRODir0[,case]
    tmpRO0 = knownRO0[case]
    tmpPO0 = knownPOminRO0[case]
    init.min = min(which(tmpdir0==1))
    init.max = max(which(tmpdir0==1))

    if (init.min>1) {
      for (i in (init.min-1):1) {
        if (plat.names2[i] %in% rownames(intron.nonout)[which(intron.nonout[,case]==0)]) {
          next
        }
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,qrsc=TRUE)[case]
        tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }

    if (init.max<length(tmpdir0)) {
      for (i in (init.max+1):length(tmpdir0)) {
        if (plat.names2[i] %in% rownames(intron.nonout)[which(intron.nonout[,case]==0)]) {
          next
        }
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,qrsc=TRUE)[case]
        tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }
    if (sum((knownRODir0[,case]-tmpdir0)^2)>0) {
      knownRODir[,case] = tmpdir0
      # colnames(knownRODir)[case] = strsplit(rownames(knownRODir)[min(which(knownRODir[,case]==1))],"[.]")[[1]][1]
    }
    rm(tmpdir0)
    cat(case,"\n")
  }
  knownPOminRO = diag(get_POgivenB(X=normProjData,B=size.adjMat%*%knownRODir,qrsc=TRUE))
  knownRO = diag(get_resdZsum(X=normProjData,B=size.adjMat%*%knownRODir,L1=FALSE)) # Residual (Z-score) outlyingness
  # knownPOminRO = diag(knownPOmat[apply(knownROmat,2,which.min),])


  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(knownPOminRO^2 < 1e-10)
  for (df in 1:20) {
    # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
    temp.out=which(abs(knownPOminRO)>sqrt(qchisq(0.001,df=df,lower.tail=F)))
    ks.output=ks.test(unique((knownPOminRO^2)[which(! c(1:length(knownPOminRO)) %in% c(temp.out,temp.zero))]),pchisq,df)
    ks.pval[df]=ks.output$p.value
    ks.stat[df]=ks.output$statistic
    # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
    #           round(ks.stat[df],digits=3)),"\n")
  }
  df=which.min(ks.stat)
  # cutoff = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),6)
  cutoff = sqrt(qchisq(p=(1-siglev),df=df))

  # Save results
  OStmp = knownPOminRO
  MODtmp = knownRODir
  OUTtmp = which(abs(OStmp)>cutoff)

  ## 2. Detecting outliers from ATS/ATT directions
  atDir = build_atDir(plat.table=plat.table,
                      inputData=normData,Ranges=Ranges,JSR.table=JSR.table,
                      cutoff=cutoff)
  if (ncol(atDir)>0) {
    adj.atDir = size.adjMat%*%atDir
    # atPOmat = get_POgivenB(X=normProjData,B=adj.atDir,qrsc=TRUE)
    # atPO = diag(atPOmat[apply(abs(atPOmat),2,which.max),])
    # atPODir = atDir[,apply(atPOmat,2,which.max)]
    atROmat = get_resdZsum(X=normProjData,B=adj.atDir,L1=FALSE,qrsc=TRUE)
    atROmat0 = atROmat
    for (ii in 1:nrow(atROmat)) {
      tmp.name = gsub("ATS","I",unlist(strsplit(rownames(atROmat)[ii],"[.]"))[1])
      tmp.name = gsub("ATT","I",tmp.name)
      atROmat0[ii,which(intron.nonout[which(rownames(intron.nonout)==tmp.name),]==0)] = 1000
    }
    atRO = diag(atROmat[apply(atROmat0,2,which.min),])
    atRODir = atDir[,apply(atROmat0,2,which.min)]
    atPOminRO = diag(get_POgivenB(X=normProjData,B=size.adjMat%*%atRODir,qrsc=TRUE))
    # atPOminRO = diag(atPOmat[apply(atROmat0,2,which.min),])
    # atSC = which((abs(atPOminRO)>cutoff) & ((abs(atPOminRO)>abs(OStmp)) & (atRO<knownRO)))
    atSC = which((atPOminRO>cutoff) & (atRO<knownRO))
    atSC = atSC[which((knownRO^2 - atRO^2)[atSC]>qchisq(p=0.01,df=1,lower.tail=F))]

    ## Update the OS and MOD
    OStmp[atSC] = atPOminRO[atSC]
    MODtmp[,atSC] = atRODir[,atSC]
    colnames(MODtmp)[atSC] = colnames(atRODir)[atSC]
  }


  ## Get final results table
  # Find region tags
  cite_first = function(t) {
    return(strsplit(t,"[.]")[[1]][1])
  }
  find_region.tag = function(MOD) {
    region.names = rownames(MOD)
    inner_fn = function(x) {
      paste(unique(sapply(region.names[which(! x==0)],cite_first)),collapse=",")
    }
    # region.tag1 = sapply(colnames(MOD),cite_first)
    # region.tag2 = apply(MOD,2,inner_fn)
    # region.tag0 = paste(region.tag1,region.tag2,sep=",")
    # region.tag = sapply(region.tag0,FUN=function(t){paste(unique(unlist(strsplit(t,","))),collapse=",")})
    region.tag = apply(MOD,2,inner_fn)
    return(region.tag)
  }
  region.tag = find_region.tag(MOD=MODtmp)
  # paste(unique(sapply(region.names[which(! MODtmp[,4]==0)],cite_first)),collapse=",")

  # Results
  outliers = which(abs(OStmp)>cutoff)
  output.tag = region.tag[outliers]
  basis.tag = sapply(colnames(MODtmp)[outliers],cite_first)
  # output.tag = colnames(MODtmp[,outliers])
  # if (length(which(grepl("[.]",output.tag)))>0) {
  #   output.tag[which(grepl("[.]",output.tag))] =
  #     sapply(output.tag[which(grepl("[.]",output.tag))],FUN=function(t){strsplit(t,"[.]")[[1]][1]})
  # }
  globalResult.table = data.frame(Outlier=outliers,
                                  Statistic=round(abs(OStmp[outliers]),digits=3),
                                  Region=locate_region(tag=colnames(MODtmp[,outliers]),Ranges=Ranges,JSR.table=JSR.table),
                                  Region.tag=output.tag,
                                  Basis.tag=basis.tag)

  ## Output
  outputObject = list(table=globalResult.table,
                      Outlier=globalResult.table$Outlier,
                      OS=abs(OStmp),
                      signOS=OStmp,
                      platBasis=plat.baseMat,
                      platMOD=MODtmp,
                      cutoff=cutoff)
  class(outputObject) = append(class(outputObject),"GSCOutput")
  return(outputObject)
}

miscGlobal_test2 = function(inputData,Ranges,JSR.table,
                            siglev=1e-04,ADcutoff=3) {

  n = dim(inputData)[2]
  plat.junctions = plat_gene(Ranges=Ranges,JSR.table=JSR.table)
  plat.table = build_platTable(Ranges=Ranges,JSR.table=JSR.table)
  plat.baseMat = build_baseMat(plat.table)
  size.adjMat = sqrt(t(plat.baseMat)%*%plat.baseMat)

  ## 1. Detecting outliers from known directions
  knownDir = build_knownDir(plat.table=plat.table,Ranges=Ranges,JSR.table=JSR.table)
  platBasisDir = apply(plat.baseMat,2,get_unitdir)
  KnownBasisDir = apply(plat.baseMat%*%knownDir,2,get_unitdir)

  normProjData = t(platBasisDir) %*% normData # low-dimensional data object to be used for outlier detection
  adj.knownDir = size.adjMat%*%knownDir # adjust size of plats

  # Get primary projection outlyingness
  knownPOmat = get_POgivenB(X=normProjData,B=adj.knownDir,qrsc=TRUE)
  knownPO = diag(knownPOmat[apply(abs(knownPOmat),2,which.max),])
  # knownPODir = knownDir[,apply(abs(knownPOmat),2,which.max)]

  # Get sum of squares of the residuals
  knownROmat = get_resdZsum(X=normProjData,B=adj.knownDir,L1=FALSE) # Residual (Z-score) outlyingness
  knownRO0 = diag(knownROmat[apply(knownROmat,2,which.min),])
  knownRODir0 = knownDir[,apply(knownROmat,2,which.min)]
  # knownPOminRO0 = diag(knownPOmat[apply(knownROmat,2,which.min),])
  knownPOminRO0 = diag(get_POgivenB(X=normProjData,B=size.adjMat%*%knownRODir0,qrsc=TRUE))
  knownRODir = knownRODir0

  for (case in 1:n) {
    tmpdir0 = knownRODir0[,case]
    tmpRO0 = knownRO0[case]
    tmpPO0 = knownPOminRO0[case]
    init.min = min(which(tmpdir0==1))
    init.max = max(which(tmpdir0==1))

    if (init.min>1) {
      for (i in (init.min-1):1) {
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,qrsc=TRUE)[case]
        tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }

    if (init.max<length(tmpdir0)) {
      for (i in (init.max+1):length(tmpdir0)) {
        tmpdir1 = tmpdir0
        tmpdir1[i] = 1
        direction1 = get_unitdir(size.adjMat%*%tmpdir1)
        tmpRO1 = get_resdZsum1d(X=normProjData,direction=direction1,qrsc=TRUE)[case]
        tmpPO1 = get_POgivenB(X=normProjData,B=direction1,qrsc=TRUE)[case]
        if ((tmpRO0^2 - tmpRO1^2)>qchisq(p=0.001,df=1,lower.tail=F)) {
          tmpdir0 = tmpdir1
          tmpRO0 = tmpRO1
          tmpPO0 = tmpPO1
          # cat(paste(i,"added"),"\n")
        }
      }
    }
    if (sum((knownRODir0[,case]-tmpdir0)^2)>0) {
      knownRODir[,case] = tmpdir0
      # colnames(knownRODir)[case] = strsplit(rownames(knownRODir)[min(which(knownRODir[,case]==1))],"[.]")[[1]][1]
    }
    rm(tmpdir0)
    cat(case,"\n")
  }
  knownPOminRO = diag(get_POgivenB(X=normProjData,B=size.adjMat%*%knownRODir,qrsc=TRUE))
  knownRO = diag(get_resdZsum(X=normProjData,B=size.adjMat%*%knownRODir,L1=FALSE)) # Residual (Z-score) outlyingness
  # knownPOminRO = diag(knownPOmat[apply(knownROmat,2,which.min),])

  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(knownPOminRO^2 < 1e-10)
  for (df in 1:20) {
    # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
    temp.out=which(abs(knownPOminRO)>sqrt(qchisq(0.001,df=df,lower.tail=F)))
    ks.output=ks.test((knownPOminRO^2)[which(! c(1:length(knownPOminRO)) %in% c(temp.out,temp.zero))],pchisq,df)
    ks.pval[df]=ks.output$p.value
    ks.stat[df]=ks.output$statistic
    # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
    #           round(ks.stat[df],digits=3)),"\n")
  }
  df=which.min(ks.stat)
  cutoff = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),5)
  # cutoff = sqrt(qchisq(p=(1-siglev),df=11))

  # Save results
  OStmp = knownPOminRO
  MODtmp = knownRODir
  OUTtmp = which(abs(OStmp)>cutoff)

  ## 2. Detecting outliers from ATS/ATT directions
  atDir = build_atDir(plat.table=plat.table,
                      inputData=normData,Ranges=Ranges,
                      JSR.table=JSR.table,
                      cutoff=cutoff)
  adj.atDir = sqrt(t(plat.baseMat)%*%plat.baseMat)%*%atDir
  atPOmat = get_POgivenB(X=normProjData,B=adj.atDir,qrsc=TRUE)
  # atPO = diag(atPOmat[apply(abs(atPOmat),2,which.max),])
  atPODir = atDir[,apply(abs(atPOmat),2,which.max)]
  atROmat = get_resdZsum(X=normProjData,B=adj.atDir,L1=FALSE,qrsc=TRUE)
  atRO = diag(atROmat[apply(atROmat,2,which.min),])
  atPOminRO = diag(atPOmat[apply(atROmat,2,which.min),])
  # atSC = which((abs(atPOminRO)>cutoff) & ((abs(atPOminRO)>abs(OStmp)) & (atRO<knownRO)))
  atSC = which((abs(atPOminRO)>cutoff) & (atRO<knownRO))
  atSC = atSC[which((knownRO^2 - atRO^2)[atSC]>qchisq(p=0.01,df=1,lower.tail=F))]

  ## Update the OS and MOD
  OStmp[atSC] = atPOminRO[atSC]
  MODtmp[,atSC] = atPODir[,atSC]
  colnames(MODtmp)[atSC] = colnames(atPODir)[atSC]

  ## Get final results table
  # Find region tags
  cite_first = function(t) {
    return(strsplit(t,"[.]")[[1]][1])
  }
  find_region.tag = function(MOD) {
    region.names = rownames(MOD)
    inner_fn = function(x) {
      paste(unique(sapply(region.names[which(! x==0)],cite_first)),collapse=",")
    }
    # region.tag1 = sapply(colnames(MOD),cite_first)
    # region.tag2 = apply(MOD,2,inner_fn)
    # region.tag0 = paste(region.tag1,region.tag2,sep=",")
    # region.tag = sapply(region.tag0,FUN=function(t){paste(unique(unlist(strsplit(t,","))),collapse=",")})
    region.tag = apply(MOD,2,inner_fn)
    return(region.tag)
  }
  region.tag = find_region.tag(MOD=MODtmp)
  # paste(unique(sapply(region.names[which(! MODtmp[,4]==0)],cite_first)),collapse=",")

  # Results
  outliers = which(abs(OStmp)>cutoff)
  output.tag = region.tag[outliers]
  basis.tag = sapply(colnames(MODtmp)[outliers],cite_first)
  # output.tag = colnames(MODtmp[,outliers])
  # if (length(which(grepl("[.]",output.tag)))>0) {
  #   output.tag[which(grepl("[.]",output.tag))] =
  #     sapply(output.tag[which(grepl("[.]",output.tag))],FUN=function(t){strsplit(t,"[.]")[[1]][1]})
  # }
  globalResult.table = data.frame(Outlier=outliers,
                                  Statistic=round(abs(OStmp[outliers]),digits=3),
                                  Region=locate_region(tag=colnames(MODtmp[,outliers]),Ranges=Ranges,JSR.table=JSR.table),
                                  Region.tag=output.tag,
                                  Basis.tag=basis.tag)

  ## Output
  outputObject = list(table=globalResult.table,
                      Outlier=globalResult.table$Outlier,
                      OS=abs(OStmp),
                      signOS=OStmp,
                      platBasis=plat.baseMat,
                      platMOD=MODtmp,
                      cutoff=cutoff)
  class(outputObject) = append(class(outputObject),"GSCOutput")
  return(outputObject)
}


miscLocal_test2 = function(miscGlobalResult,
                           inputData,pileupData,
                           Ranges,JSR.table,
                           siglev=1e-4,cutoff=NULL,
                           ADcutoff=4,windowSize=50) {
  n = ncol(inputData); d = nrow(inputData);
  if (is.null(cutoff)) {
    cutoff.g = miscGlobalResult$cutoff
  }
  out1 = miscGlobalResult$Outlier
  curIDs = which(! c(1:n) %in% out1)
  inputData2 = inputData[,curIDs]
  pileupData2 = pileupData[,curIDs]
  n2 = dim(pileupData2)[2]
  
  get_rps = function(X,dir) {
    # robust projection scores by projecting X onto dir.
    require(MASS)
    return(apply(X,2,FUN=function(t){rlm(t ~ dir - 1, psi=psi.bisquare)$coefficients}))
  }
  z.pca = pca.hy(inputData2,subt.mean=F)
  maincomp = function(x) {
    # 1 = main / 0 = outlier
    require(diptest)
    if (dip.test(x)$p.value<0.005) {
      (max(abs(correct_PO(pd.rate.hy(x,qrsc=T))))<3)
    } else {
      (max(abs(correct_PO(pd.rate.hy(x,qrsc=T))))<4)
    }
  }
  pc.j = which(apply(z.pca$projmat[1:10,],1,maincomp)==1)
  if (length(pc.j)>0) {
    rprojmat = apply(matrix(z.pca$dirmat[,pc.j],ncol=length(pc.j)),2,
                     FUN=function(t){get_rps(X=inputData2,dir=t)})
    inputData2 = inputData2 - z.pca$dirmat[,pc.j]%*%t(rprojmat)
  } else {
    inputData2 = inputData2
  }
  
  ## 1. Detect outliers from cryptic junction directions
  crypBasisRes = build_crypticDir(Ranges=Ranges,JSR.table=JSR.table)
  crypBasisDir = crypBasisRes$BasisDir
  crypJDir.names = colnames(crypBasisDir)
  crypJDir.tag = sapply(crypJDir.names,
                      FUN=function(t) {as.character(JSR.table$Region.tag)[which(as.character(JSR.table$JV.tag)==t)]})

  crypOS = rep(0,n)
  if (length(crypJDir.names)>0) {
    ## 1-1. Get projection outlyingness
    crypPOmat = matrix(get_POgivenB2(X=inputData2,B=crypBasisDir,qrsc=TRUE,ADcutoff=ADcutoff),ncol=n2)

    ## 1-2. Non-unimodality / Extreme skewness
    crypPOmat = sweep(x=crypPOmat,1,STATS=apply(crypPOmat,1,filter_multimodal),"*")

    ## 1-3. Non-outliers from intronic regions
    intron.nonout2 = matrix(1,nrow=ncol(crypBasisDir),ncol=n2)
    for (i in which(grepl(pattern="I",crypJDir.tag)==T)) {
      carea = which(crypBasisDir[,i]^2>0)
      ie = as.numeric(unlist(strsplit(crypJDir.tag[i],"I"))[2])
      carea1 = c(Ranges$lRanges[ie,2]:Ranges$lRanges[ie,3])
      carea2 = c(Ranges$lRanges[(ie+1),2]:Ranges$lRanges[(ie+1),3])
      intron.nonout2[i,] = apply(pileupData2,2,
                                 FUN=function(t){filter_intron_onoff(x=t,carea=carea,carea1=carea1,carea2=carea2,cryptic=TRUE)})
    }

    ## 1-4. Outlier detection
    crypPO = diag(crypPOmat[apply(abs(crypPOmat)*intron.nonout2,2,which.max),])
    crypPODir0 = crypBasisDir[,apply(abs(crypPOmat)*intron.nonout2,2,which.max)]
    # rows.incl = c(1:nrow(crypPOmat))[(! 1:nrow(crypPOmat) %in% which((grepl("E",get_Tag(JSR.table$LBE.position[sapply(rownames(crypPOmat),FUN=function(t){which(JSR.table$JV.tag==t)})]))) & (apply(crypPOmat,1,ADstatWins.hy)>ADcutoff)))]
    # crypPOmat = crypPOmat[rows.incl,]
    # crypPO = diag(crypPOmat[apply(abs(crypPOmat),2,which.max),])
    # crypPODir0 = crypBasisDir[,rows.incl[apply(abs(crypPOmat),2,which.max)]]

    crypOS0 = abs(crypPO)
    ks.pval=ks.stat=rep(0,20)
    temp.zero=which(crypOS0^2 < 1e-10)
    if (length(temp.zero)<(0.5*length(crypOS0))) {
      for (df in 1:20) {
        # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
        temp.out=which(crypOS0>sqrt(qchisq(0.001,df=df,lower.tail=F)))
        ks.output=ks.test(unique(crypOS0[which(! c(1:length(crypOS0)) %in% c(temp.out,temp.zero))])^2,pchisq,df)
        ks.pval[df]=ks.output$p.value
        ks.stat[df]=ks.output$statistic
        # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
        #           round(ks.stat[df],digits=3)),"\n")
      }
      df=which.min(ks.stat)
      # cutoff = sqrt(qchisq(siglev,df=df,lower.tail=F))
      cutoff = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),4)
    } else {
      cutoff = cutoff.g
    }
    crypOS[curIDs] = crypPO
  } else {
    crypPODir0 = matrix(0,ncol=n2,nrow=nrow(crypBasisDir))
    colnames(crypPODir0) = rep(".",n2)
  }

  crypPODir = matrix(0,ncol=n,nrow=nrow(crypPODir0))
  colnames(crypPODir) = rep(".",n)
  colnames(crypPODir)[which(! c(1:n) %in% out1)] = colnames(crypPODir0)
  crypPODir[,which(! c(1:n) %in% out1)] = crypPODir0
  out2 = which(abs(crypOS)>cutoff)

  basis.tag = colnames(crypPODir)[out2]
  region.tag = sapply(colnames(crypPODir)[out2],
                      FUN=function(t) {as.character(JSR.table$Region.tag)[which(as.character(JSR.table$JV.tag)==t)]})
  crypResult = data.frame(Outlier=out2,
                          Statistic=round(abs(crypOS[out2]),digits=3),
                          Region=sapply(colnames(crypPODir)[out2],FUN=function(t){crypBasisRes$BasisDir.pos[which(crypJDir.names==t)]}),
                          Region.tag=region.tag,
                          Basis.tag=basis.tag)

  ## 2. Detect outliers from window directions
  out1 = c(miscGlobalResult$Outlier,crypResult$Outlier)
  exonset = Ranges$lRanges
  if (nrow(exonset)==1) {
    exon.base=c(exonset[1,2]:exonset[1,3])
  } else {
    exon.base=c()
    for (i in 1:(nrow(exonset)-1)) {
      exon.base=c(exon.base,c(exonset[i,2]:exonset[i,3]))
    }
  }

  residualData2 = inputData2[,which(! colnames(inputData2) %in% colnames(inputData)[out1])]
  pileupData2 = pileupData[,which(! colnames(pileupData) %in% colnames(inputData)[out1])]
  
  onoff_res = get_offstat(residualData=residualData2,pileupData=pileupData2,
                          exonset=exonset,ADcutoff=ADcutoff,
                          windowSize=windowSize,readconstr=10)
  # localOS0 = onoff_res$stat
  onoff_res$tMOD = -onoff_res$tMOD
  tmpIDs = which(onoff_res$stat==0)
  localPOmat = matrix(0,ncol=ncol(residualData2),nrow=ncol(residualData2))
  localPOmat = get_POgivenB2(X=residualData2,B=onoff_res$tMOD,qrsc=TRUE,ADcutoff=ADcutoff)
  localPOmat = sweep(x=localPOmat,1,STATS=apply(localPOmat,1,filter_multimodal),"*")
  localPOmat[tmpIDs,tmpIDs] = 0

  localOS0 = rep(0,n)
  localMOD0 = matrix(0,nrow=d,ncol=n)
  localOS0[which(! c(1:n) %in% out1)] = diag(localPOmat)
  localMOD0[,which(! c(1:n) %in% out1)] = onoff_res$tMOD

  localOS = crypOS
  localOS[which(abs(localOS0)>abs(crypOS))] = localOS0[which(abs(localOS0)>abs(crypOS))]
  localOS[out2] = crypOS[out2]
  localMOD = crypPODir
  localMOD[,which(abs(localOS0)>abs(crypOS))] = localMOD0[,which(abs(localOS0)>abs(crypOS))]
  colnames(localMOD)[which(abs(localOS0)>abs(crypOS))] = "."

  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(localOS^2 < 1e-10)
  if (length(temp.zero)<(0.5*length(localOS))) {
    for (df in 1:20) {
      # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
      temp.out=which(localOS^2>qchisq(0.001,df=df,lower.tail=F))
      ks.output=ks.test(unique(localOS[which(! c(1:length(localOS0)) %in% c(temp.out,temp.zero))])^2,pchisq,df)
      ks.pval[df]=ks.output$p.value
      ks.stat[df]=ks.output$statistic
      # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
      #           round(ks.stat[df],digits=3)),"\n")
    }
    df=which.min(ks.stat)
    cutoff = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),cutoff)
  }
  out3 = which(abs(localOS0)>cutoff)

  if (length(out3)>0) {
    localOutRegion = apply(matrix(localMOD[,out3],ncol=length(out3)),2,FUN=function(t){collapse_junction(c(min(which(t^2>0)),max(which(t^2>0))))})
    region.tag = paste("E",find_exon(localOutRegion,Ranges=Ranges),sep="")
    basis.tag = rep(".",length(out3))
  } else {
    localOutRegion = region.tag = basis.tag = NULL
  }

  localResult.table = data.frame(Outlier=out3,
                                 Statistic=round(abs(localOS[out3]),digits=3),
                                 Region=localOutRegion,
                                 Region.tag=region.tag,
                                 Basis.tag=basis.tag)
  localResult.table = rbind(localResult.table,
                            crypResult[which(crypResult$Statistic>cutoff),])
  localResult.table = localResult.table[order(localResult.table$Outlier),]
  rownames(localResult.table) = 1:nrow(localResult.table)
  outputObject = list(table=localResult.table,
                      Outlier=localResult.table$Outlier,
                      OS=abs(localOS),
                      signOS=localOS,
                      MOD=localMOD,
                      cutoff=cutoff)
  class(outputObject) = append(class(outputObject),"LSCOutput")
  return(outputObject)
}


# get_intron_nonout = function(inputData,pileupData,Ranges) {
#   n = ncol(inputData)
#   exonset = Ranges$lRanges
#   nexon = nrow(exonset)
#
#   intron_onstat_inner=function(x,carea,carea1,carea2){
#     readconstr=10
#     lcarea=length(carea)
#     y=max(quantile(x[carea[6:ceiling(lcarea/3)]],probs=0.75),
#           quantile(x[carea[ceiling(lcarea/3):(2*ceiling(lcarea/3))]],probs=0.75),
#           quantile(x[carea[(2*ceiling(lcarea/3)):(lcarea-6)]],probs=0.75))
#     max_exonread = max(quantile(x[carea1],probs=0.5),quantile(x[carea2],probs=0.5));
#     if (max_exonread>readconstr) {
#       readconstr2=max(5,0.1*max_exonread)
#     } else {
#       readconstr2=readconstr
#     }
#     if (max_exonread>100) {
#       ((y>=readconstr) & (y>=readconstr2))
#     } else {
#       ((y>=readconstr) | (y>=readconstr2))
#     }
#   }
#
#   medianmat = matrix(0,nrow=(nexon-1),ncol=n)
#   if (nexon > 1) {
#     for (i in 1:(nexon-1)) {
#       carea = c((exonset[i,3]+1):(exonset[(i+1),2]-1))
#       if (length(carea)>20) {
#         carea1 = c(exonset[i,2]:exonset[i,3])
#         carea2 = c(exonset[(i+1),2]:exonset[(i+1),3])
#         medianmat[i,] = apply(pileupData,2,
#                               FUN=function(t){intron_onstat_inner(x=t,carea=carea,carea1=carea1,carea2=carea2)})
#       }
#     }
#   }
#   rownames(medianmat) = paste("I",1:dim(medianmat)[1],sep="")
#   colnames(medianmat) = colnames(inputData)
#   return(medianmat)
# }

# get_intron_nonout = function(pileupData,Ranges) {
#   n = ncol(pileupData)
#   exonset = Ranges$lRanges
#   nexon = nrow(exonset)
#
#   intron_onstat_inner=function(x,carea,carea1,carea2){
#     lcarea=length(carea)
#     max_exonread = max(quantile(x[carea1[c(floor(0.5*length(carea1)):length(carea1))]],probs=0.5),
#                        quantile(x[carea2[c(1:floor(0.5*length(carea2)))]],probs=0.5))
#     if (max_exonread<5) {
#       return(0)
#     } else {
#       y=max(quantile(x[carea[6:ceiling(lcarea/3)]],probs=0.75),
#             quantile(x[carea[ceiling(lcarea/3):(2*ceiling(lcarea/3))]],probs=0.75),
#             quantile(x[carea[(2*ceiling(lcarea/3)):(lcarea-6)]],probs=0.75))
#       cutoff_ratio = IR_cutoff_fn(max_exonread)
#       return((y>=cutoff_ratio*max_exonread))
#     }
#   }
#   medianmat = matrix(0,nrow=(nexon-1),ncol=n)
#   if (nexon > 1) {
#     for (i in 1:(nexon-1)) {
#       carea = c((exonset[i,3]+1):(exonset[(i+1),2]-1))
#       if (length(carea)>20) {
#         carea1 = c(exonset[i,2]:exonset[i,3])
#         carea2 = c(exonset[(i+1),2]:exonset[(i+1),3])
#         medianmat[i,] = apply(pileupData,2,
#                               FUN=function(t){intron_onstat_inner(x=t,carea=carea,carea1=carea1,carea2=carea2)})
#       }
#     }
#   }
#   rownames(medianmat) = paste("I",1:dim(medianmat)[1],sep="")
#   colnames(medianmat) = colnames(pileupData)
#   return(medianmat)
# }

# get_intron_nonout = function(inputData,pileupData,Ranges) {
#   n = ncol(inputData)
#   exonset = Ranges$lRanges
#   nexon = nrow(exonset)
#
#   intron_onstat_inner=function(x,carea,carea1,carea2){
#     readconstr = 10
#     lcarea=length(carea)
#     y=max(quantile(x[carea[6:ceiling(lcarea/3)]],probs=0.75),
#           quantile(x[carea[ceiling(lcarea/3):(2*ceiling(lcarea/3))]],probs=0.75),
#           quantile(x[carea[(2*ceiling(lcarea/3)):(lcarea-6)]],probs=0.75))
#     max_exonread = max(quantile(x[carea1],probs=0.5),quantile(x[carea2],probs=0.5))
#
#     if (max_exonread>readconstr) {
#       readconstr2=max(5,0.1*max_exonread)
#     } else {
#       readconstr2=readconstr
#     }
#     if (max_exonread>100) {
#       ((y>=readconstr) & (y>=readconstr2))
#     } else {
#       ((y>=readconstr) | (y>=readconstr2))
#     }
#   }
#
#   medianmat = matrix(0,nrow=(nexon-1),ncol=n)
#   if (nexon > 1) {
#     for (i in 1:(nexon-1)) {
#       carea = c((exonset[i,3]+1):(exonset[(i+1),2]-1))
#       if (length(carea)>20) {
#         carea1 = c(exonset[i,2]:exonset[i,3])
#         carea2 = c(exonset[(i+1),2]:exonset[(i+1),3])
#         medianmat[i,] = apply(pileupData,2,
#                               FUN=function(t){intron_onstat_inner(x=t,carea=carea,carea1=carea1,carea2=carea2)})
#       }
#     }
#   }
#   rownames(medianmat) = paste("I",1:dim(medianmat)[1],sep="")
#   colnames(medianmat) = colnames(inputData)
#   return(medianmat)
# }
