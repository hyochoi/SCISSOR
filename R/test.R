#' Get final results table from SCISSOR
#'
#' @export
get_finalTable = function(miscGlobalResult,miscLocalResult,jsResult) {

  region.tag = as.character(jsResult$Region.tag)
  region.tag[which(region.tag==".")] = as.character(jsResult$JV.tag)[which(region.tag==".")]
  jsResult$Region.tag = region.tag

  miscGlobalResult$table$Type = rep("Global",dim(miscGlobalResult$table)[1])
  miscLocalResult$table$Type = rep("Local",dim(miscLocalResult$table)[1])
  Stable = rbind(miscGlobalResult$table,miscLocalResult$table) # Scissor results table
  Jtable = jsResult %>% select(Outlier, VJF, JV.class, Region.tag, JV.tag) # Junction results table
  merged0 = merge(x=Jtable,y=Stable,by=c("Outlier"),all.x=T)
  # Jtable = merged0[which(apply(merged0 %>% select(Region.tag.x,Region.tag.y),1,
  #                              FUN=function(t){grepl(pattern=t[1],t[2])})==1),] %>% select(Outlier,VJF,JV.class,Region.tag.x,JV.tag)
  Jtable = merged0[which((apply(merged0 %>% select(Region.tag.x,Region.tag.y),1,
                                FUN=function(t){grepl(pattern=t[1],t[2])})==1) | (grepl("J",as.character(merged0$Region.tag.x)))),] %>% select(Outlier,VJF,JV.class,Region.tag.x,JV.tag)
  colnames(Jtable)[which(colnames(Jtable)=="Region.tag.x")] = "Region.tag"

  # merged = merge(x=Stable,y=Jtable,by=c("Outlier","Region.tag"),all.x=T)
  merged = merge(x=Stable,y=Jtable,by=c("Outlier"),all.x=T)
  colnames(merged)[which(colnames(merged)=="Region.tag.x")] = "Region.tag"
  colnames(merged)[which(colnames(merged)=="Region.tag.y")] = "JV.Region.tag"

  # merged[which((! apply(merged %>% select(JV.Region.tag,Region.tag),1,
  #                FUN=function(t){grepl(pattern=t[1],t[2])})==1) & (grepl("J",as.character(merged$JV.Region.tag)))),]
  ## Classify outlier types
  class.g = as.character(merged$JV.class)
  for (i in which(class.g=="TBD")) {
    case = merged$Outlier[i]
    tmp.tag = as.character(merged$Region.tag)[i]
    if (grepl("E",tmp.tag)) {
      if (miscGlobalResult$signOS[case]<0) {
        class.g[i] = "ES"
      } else if (miscGlobalResult$signOS[case]>0) {
        class.g[i] = "AE"
      }
    } else if (grepl("J",tmp.tag)) {
      if (miscGlobalResult$signOS[case]<0) {
        class.g[i] = "Del"
      } else if (miscGlobalResult$signOS[case]>0) {
        class.g[i] = "Unknown"
      }
    } else {
      class.g[i] = "Unknown"
    }
  }
  for (i in which(is.na(class.g)==T)) {
    case = merged$Outlier[i]
    tmp.tag = as.character(merged$Region.tag)[i]
    if (grepl("ATS",tmp.tag)) {
      class.g[i] = "ATS"
    } else if (grepl("ATT",tmp.tag)) {
      class.g[i] = "ATT"
    } else if (grepl("I",tmp.tag)) {
      class.g[i] = "IR"
    } else if (grepl("E",tmp.tag)) {
      if (as.character(merged$Type)[i]=="Local") {
        class.g[i] = "Del"
      } else {
        if (miscGlobalResult$signOS[case]<0) {
          class.g[i] = "Del"
        } else if (miscGlobalResult$signOS[case]>0) {
          class.g[i] = "Unknown"
        }
      }
    } else {
      class.g[i] = "Unknown"
    }
  }
  merged$Class = class.g

  ## VJF for ATT events
  class.tag = as.character(merged$Class)
  for (i in which((class.tag=="ATS") | (class.tag=="ATT"))) {
    merged[i,]
    case = merged$Outlier[i]
    tmp.tag = paste("I",strsplit(as.character(merged$Region.tag)[i],class.tag[i])[[1]][2],sep="")
    merged[i,c(6,8)] = Jtable[which((Jtable$Outlier==case) & (as.character(Jtable$Region.tag)==tmp.tag)),] %>% select(VJF, JV.tag)
  }

  merged = merged %>% select(Outlier, Statistic, Class, Region, Region.tag, VJF, JV.Region.tag, JV.tag, Type)
  return(merged)
}

#' Detect local shape changes
#'
#' @export
miscGlobal_test = function(inputData,Ranges,JSR.table,
                           siglev=1e-04,ADcutoff=3) {

  plat.junctions = plat_gene(Ranges=Ranges,JSR.table=JSR.table)
  plat.table = build_platTable(Ranges=Ranges,JSR.table=JSR.table)
  plat.baseMat = build_baseMat(plat.table)

  ## 1. Detecting outliers from known directions
  knownDir = build_knownDir(plat.table=plat.table,Ranges=Ranges,JSR.table=JSR.table)
  platBasisDir = apply(plat.baseMat,2,get_unitdir)
  KnownBasisDir = apply(plat.baseMat%*%knownDir,2,get_unitdir)

  normProjData = t(platBasisDir) %*% normData # low-dimensional data object to be used for outlier detection
  adj.knownDir = sqrt(t(plat.baseMat)%*%plat.baseMat)%*%knownDir # adjust size of plats

  # Get primary projection outlyingness
  knownPOmat = get_POgivenB(X=normProjData,B=adj.knownDir)
  knownPO = diag(knownPOmat[apply(abs(knownPOmat),2,which.max),])
  # knownPODir = knownDir[,apply(abs(knownPOmat),2,which.max)]

  # Get sum of squares of the residuals
  knownROmat = get_resdZsum(X=normProjData,B=adj.knownDir,L1=FALSE) # Residual (Z-score) outlyingness
  knownRO = diag(knownROmat[apply(knownROmat,2,which.min),])
  knownRODir = knownDir[,apply(knownROmat,2,which.min)]
  knownPOminRO = diag(knownPOmat[apply(knownROmat,2,which.min),])

  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(knownPO^2 < 1e-10)
  for (df in 1:20) {
    # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
    temp.out=which(knownPO>sqrt(qchisq(0.001,df=df,lower.tail=F)))
    ks.output=ks.test(knownPO[which(! c(1:length(knownPO)) %in% c(temp.out,temp.zero))]^2,pchisq,df)
    ks.pval[df]=ks.output$p.value
    ks.stat[df]=ks.output$statistic
    # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
    #           round(ks.stat[df],digits=3)),"\n")
  }
  df=which.min(ks.stat)
  cutoff = sqrt(qchisq(siglev,df=df,lower.tail=F))
  # cutoff = sqrt(qchisq(p=(1-siglev),df=11))

  # Save results
  OStmp = knownPOminRO
  MODtmp = knownRODir
  OUTtmp = which(abs(OStmp)>cutoff)

  ## 2. Detecting outliers from ATS/ATT directions
  atDir = build_atDir(plat.table=plat.table,inputData=normData,Ranges=Ranges,JSR.table=JSR.table,
                      cutoff=cutoff)
  adj.atDir = sqrt(t(plat.baseMat)%*%plat.baseMat)%*%atDir
  atPOmat = get_POgivenB(X=normProjData,B=adj.atDir)
  # atPO = diag(atPOmat[apply(abs(atPOmat),2,which.max),])
  atPODir = atDir[,apply(abs(atPOmat),2,which.max)]
  atROmat = get_resdZsum(X=normProjData,B=adj.atDir,L1=FALSE)
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
  outliers = which(abs(OStmp)>cutoff)
  output.tag = colnames(MODtmp[,outliers])
  if (length(which(grepl("[.]",output.tag)))>0) {
    output.tag[which(grepl("[.]",output.tag))] =
      sapply(output.tag[which(grepl("[.]",output.tag))],FUN=function(t){strsplit(t,"[.]")[[1]][1]})
  }
  globalResult.table = data.frame(Outlier=outliers,
                                  Statistic=round(abs(OStmp[outliers]),digits=3),
                                  Region=locate_region(tag=colnames(MODtmp[,outliers]),Ranges=Ranges,JSR.table=JSR.table),
                                  Region.tag=output.tag)

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

#' Detect local shape changes
#'
#' This function discovers outlying subjects whose RNA-seq have "local" abnormal shapes
#' and provides the most outlying window-direction for each outlier.
#'
#' @param miscGlobalResult Result from \link{miscGlobal}
#' @param inputData normalized coverage matrix from \link{process_data}
#' @param pileupData raw coverage matrix, dataI from \link{process_data}
#' @param Ranges
#' @param JSR.table
#' @param windowSize a window length. Default is 100.
#' @param siglev the significance level for a robust outlier detection. Default is 1e-5.
#' IF cutoff is specified, siglev is not used.
#' @param cutoff the cutoff value for outlying statistics.
#' If NULL, the cutoff value is computed based on the specified siglev.
#'
#' @import zoo
#' @export
miscLocal_test = function(miscGlobalResult,
                          inputData,pileupData,
                          Ranges,JSR.table,
                          siglev=1e-4,cutoff=NULL,
                          ADcutoff=3,windowSize=50) {
  n = ncol(inputData); d = nrow(inputData);
  if (is.null(cutoff)) {
    cutoff = miscGlobalResult$cutoff
  }
  ## 1. Detect outliers from cryptic junction directions
  crypBasisRes = build_crypticDir(Ranges=Ranges,JSR.table=JSR.table)
  crypBasisDir = crypBasisRes$BasisDir
  crypJDir.names = colnames(crypBasisDir)
  out1 = miscGlobalResult$Outlier
  inputData2 = inputData[,which(! c(1:n) %in% out1)]

  crypPOmat = get_POgivenB(X=inputData2,B=crypBasisDir,qrsc=TRUE)
  crypPO = diag(crypPOmat[apply(abs(crypPOmat),2,which.max),])
  crypPODir0 = crypBasisDir[,apply(abs(crypPOmat),2,which.max)]
  # rows.incl = c(1:nrow(crypPOmat))[(! 1:nrow(crypPOmat) %in% which((grepl("E",get_Tag(JSR.table$LBE.position[sapply(rownames(crypPOmat),FUN=function(t){which(JSR.table$JV.tag==t)})]))) & (apply(crypPOmat,1,ADstatWins.hy)>ADcutoff)))]
  # crypPOmat = crypPOmat[rows.incl,]
  # crypPO = diag(crypPOmat[apply(abs(crypPOmat),2,which.max),])
  # crypPODir0 = crypBasisDir[,rows.incl[apply(abs(crypPOmat),2,which.max)]]

  crypOS0 = abs(crypPO)
  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(crypOS0^2 < 1e-10)
  if (length(temp.zero)<(0.5*length(temp.zero))) {
    for (df in 1:20) {
      # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
      temp.out=which(crypOS0>sqrt(qchisq(0.001,df=df,lower.tail=F)))
      ks.output=ks.test(crypOS0[which(! c(1:length(crypOS0)) %in% c(temp.out,temp.zero))]^2,pchisq,df)
      ks.pval[df]=ks.output$p.value
      ks.stat[df]=ks.output$statistic
      # cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
      #           round(ks.stat[df],digits=3)),"\n")
    }
    df=which.min(ks.stat)
    cutoff.here = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),cutoff)
  } else {
    cutoff.here = cutoff
  }

  crypOS = rep(0,n)
  crypOS[-out1] = crypPO
  crypPODir = matrix(0,ncol=n,nrow=nrow(crypPODir0))
  colnames(crypPODir) = rep(".",n)
  colnames(crypPODir)[-out1] = colnames(crypPODir0)
  crypPODir[,-out1] = crypPODir0
  out2 = which(abs(crypOS)>cutoff.here)

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

  residualData2 = inputData[,which(! c(1:n) %in% out1)]
  pileupData2 = pileupData[,which(! c(1:n) %in% out1)]

  onoff_res = get_offstat(residualData=residualData2,pileupData=pileupData2,
                          exonset=exonset,ADcutoff=ADcutoff,
                          windowSize=windowSize,readconstr=10)

  localOS0 = onoff_res$stat
  ks.pval=ks.stat=rep(0,20)
  temp.zero=which(localOS0^2 < 1e-10)
  if (length(temp.zero)<(0.5*length(temp.zero))) {
    for (df in 1:20) {
    # hist(pchisq(onoff_stat^2,df=df,lower.tail=F),main=df)
    temp.out=which(localOS0>sqrt(qchisq(0.001,df=df,lower.tail=F)))
    ks.output=ks.test(localOS0[which(! c(1:length(localOS0)) %in% c(temp.out,temp.zero))]^2,pchisq,df)
    ks.pval[df]=ks.output$p.value
    ks.stat[df]=ks.output$statistic
    cat(paste(df,"|",round(ks.pval[df],digits=5),"|",
              round(ks.stat[df],digits=3)),"\n")
  }
    df=which.min(ks.stat)
    cutoff.here = max(sqrt(qchisq(siglev,df=df,lower.tail=F)),cutoff)
  } else {
    cutoff.here = cutoff
  }

  localOS = rep(0,n)
  localOS[-out1] = -localOS0
  localOS[out2] = crypOS[out2]
  localMOD = matrix(0,ncol=n,nrow=d)
  localMOD[,-out1] = onoff_res$tMOD
  localMOD[,out2] = crypPODir[,out2]

  out3 = which(abs(localOS)>cutoff.here)
  out3 = out3[which(! out3 %in% out2)]
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
  localResult.table = rbind(localResult.table,crypResult)
  localResult.table = localResult.table[order(localResult.table$Outlier),]

  outputObject = list(table=localResult.table,
                      Outlier=localResult.table$Outlier,
                      OS=abs(localOS),
                      signOS=localOS,
                      MOD=localMOD,
                      cutoff=cutoff.here)
  class(outputObject) = append(class(outputObject),"LSCOutput")
  return(outputObject)
}

#' Find which exon contains the position
#'
#' @export
find_exon = function(position,Ranges) {
  inner_fn = function(x) {
    tmp.pos = as.numeric(split_junction(x))
    if (length(tmp.pos)==1) {
      return(which((Ranges$lRanges[,2]<=tmp.pos[1]) & (Ranges$lRanges[,3]>=tmp.pos[1])))
    } else {
      return(which((Ranges$lRanges[,2]<=tmp.pos[1]) & (Ranges$lRanges[,3]>=tmp.pos[2])))
    }
  }
  sapply(as.character(position),inner_fn)
}

#' Collect directions for potential cryptic events
#'
#' @export
build_crypticDir = function(Ranges,JSR.table) {
  ## abnormal junction basis (do not consider cryptic events)
  nexons = dim(Ranges$lRanges)[1]
  cryptic.idx = which(grepl(pattern="Cryptic",JSR.table[,which(colnames(JSR.table)=="JV.class")])==T)
  crypJDir.names = as.character(JSR.table$JV.tag)[cryptic.idx]
  crypJDir.pos = rep("0",length(cryptic.idx))

  for (j in 1:length(cryptic.idx)) {
    i = cryptic.idx[j]
    JSR.i = JSR.table[i,]
    LBE.tab = sapply(c(split_junction(as.character(JSR.i$LBE.position))),FUN=function(t){c(strsplit(strsplit(t,":")[[1]][1],"exon")[[1]][2],strsplit(t,":")[[1]][2])})
    LBE.which = which((! LBE.tab[2,]=="out") & (! LBE.tab[2,]=="0"))
    # cat(paste(j,"|",length(LBE.which)),"\n")
    if (LBE.which==2) {
      LBE = as.numeric(LBE.tab[,LBE.which])
      if (LBE[2]<0) {
        crypJDir.pos[j] = collapse_junction(c(Ranges$lRanges[LBE[1],2],Ranges$lRanges[LBE[1],2] - LBE[2]-1))
      } else {
        crypJDir.pos[j] = collapse_junction(c(Ranges$lRanges[LBE[1],2] - LBE[2],Ranges$lRanges[LBE[1],2]))
      }
    } else {
      LBE = as.numeric(LBE.tab[,LBE.which])
      if (LBE[2]<0) {
        crypJDir.pos[j] = collapse_junction(c(Ranges$lRanges[LBE[1],3] + LBE[2]+1,Ranges$lRanges[LBE[1],3]))
      } else {
        crypJDir.pos[j] = collapse_junction(c(Ranges$lRanges[LBE[1],3],Ranges$lRanges[LBE[1],3] + LBE[2]))
      }
    }
  }

  crypBasisDir = matrix(0,nrow=d,ncol=length(crypJDir.names))
  colnames(crypBasisDir) = crypJDir.names
  for (i in 1:length(crypJDir.names)) {
    jpos = split_junction(crypJDir.pos[i])
    crypBasisDir[c(jpos[1]:jpos[2]),i] = 1
  }
  crypBasisDir = apply(crypBasisDir,2,FUN=function(t){t/sqrt(sum(t^2))})
  return(list(BasisDir=crypBasisDir,BasisDir.pos=crypJDir.pos))
}

#' Collect known directions
#'
#' @export
build_knownDir = function(plat.table,Ranges,JSR.table) {
  # plat.table = build_platTable(Ranges = Ranges,JSR.table = JSR.table)
  nexons = dim(Ranges$lRanges)[1]

  ## Type 1. sparse plat basis only for exon & intron basis
  platDir.names = sort(unique(c(sapply(rownames(plat.table),FUN=function(t){strsplit(t,split="[.]")[[1]][1]}),rownames(plat.table)[which(grepl(pattern="I",rownames(plat.table))==T)])))
  platDir = matrix(0,nrow=dim(plat.table)[1],ncol=length(platDir.names))
  rownames(platDir) = rownames(plat.table)
  colnames(platDir) = platDir.names
  for (i in 1:ncol(platDir)) {
    ij = which(rownames(plat.table)==platDir.names[i])
    if (length(ij)>0) {
      platDir[ij,i] = 1
    } else {
      platDir[which(sapply(rownames(plat.table),FUN=function(t){strsplit(t,split="[.]")[[1]][1]})==platDir.names[i]),i] = 1
    }
  }

  ## Type 2. abnormal junction basis (do not consider cryptic events)
  non.cryptic.idx = which(grepl(pattern="Cryptic",JSR.table[,which(colnames(JSR.table)=="JV.class")])==F)
  JDir = matrix(0,nrow=dim(plat.table)[1],ncol=(2*dim(JSR.table)[1]))
  rownames(JDir) = rownames(plat.table)
  JDir.names = c()
  di = 1
  for (i in non.cryptic.idx) {
    JSR.i = JSR.table[i,]
    LBE.tab = sapply(c(split_junction(as.character(JSR.i$LBE.position))),FUN=function(t){c(strsplit(strsplit(t,":")[[1]][1],"exon")[[1]][2],strsplit(t,":")[[1]][2])})
    if (sum(grepl(pattern="out",LBE.tab[2,]))>0) {
      if ((diff(as.numeric(LBE.tab[1,]))==0) & (sum(grepl(pattern="0",LBE.tab[2,]))>0)) {
        next
      }
      exon.bp.tmp = as.numeric(split_junction(as.character(JSR.i$junctions.l)))
      if (exon.bp.tmp[1]<0) {
        basei = which((as.numeric(sapply(as.character(plat.table$Range),split_junction)[2,])<exon.bp.tmp[2]) & (plat.table$Domain=="exon"))
        JDir[basei,di] = 1
        JDir.names = c(JDir.names,as.character(JSR.i$JV.tag))
        di = di + 1
        if (diff(as.numeric(LBE.tab[1,]))>(nexons/3)) {
          JDir[which(plat.table$Domain=="exon")[which(! which(plat.table$Domain=="exon") %in% basei)],di] = 1
          JDir.names = c(JDir.names,as.character(JSR.i$JV.tag))
          di = di + 1
        }
      } else {
        basei = which((as.numeric(sapply(as.character(plat.table$Range),split_junction)[1,])>exon.bp.tmp[1]) & (plat.table$Domain=="exon"))
        JDir[basei,di] = 1
        JDir.names = c(JDir.names,as.character(JSR.i$JV.tag))
        di = di + 1
        if (diff(as.numeric(LBE.tab[1,]))>(nexons/3)) {
          JDir[which(plat.table$Domain=="exon")[which(! which(plat.table$Domain=="exon") %in% basei)],di] = 1
          JDir.names = c(JDir.names,as.character(JSR.i$JV.tag))
          di = di + 1
        }
      }
    } else {
      exon.bp.tmp = as.numeric(split_junction(as.character(JSR.i$junctions.l)))
      basei = which(((as.numeric(sapply(as.character(plat.table$Range),split_junction)[2,])<exon.bp.tmp[2]) & (as.numeric(sapply(as.character(plat.table$Range),split_junction)[1,])>=exon.bp.tmp[1])) & (plat.table$Domain=="exon"))
      if (length(basei)>0) {
        JDir[basei,di] = 1
        JDir.names = c(JDir.names,as.character(JSR.i$JV.tag))
        di = di + 1
        if (diff(as.numeric(LBE.tab[1,]))>(nexons/3)) {
          JDir[which(plat.table$Domain=="exon")[which(! which(plat.table$Domain=="exon") %in% basei)],di] = 1
          JDir.names = c(JDir.names,as.character(JSR.i$JV.tag))
          di = di + 1
        }
      }
    }
    rm(basei)
  }
  JDir = JDir[,1:which.max(which(apply(JDir,2,sum)>0))]
  colnames(JDir) = JDir.names

  ## Type 3. Get SV basis
  p = dim(plat.table)[1]
  SVDir = matrix(0,nrow=p,ncol=(2*(nexons-1)))
  rownames(SVDir) = rownames(plat.table)
  colnames(SVDir) = rep(paste("SV",1:(nexons-1),sep=""),each=2)
  for (i in 1:(nexons-1)) {
    intron.j = which(sapply(rownames(plat.table),FUN=function(t){strsplit(t,"[.]")[[1]][1]})==paste0("I",i))[1]
    SVDir[which((plat.table$Domain=="exon") & (c(1:p) < intron.j)),(2*(i-1)+1)] = 1
    SVDir[which((plat.table$Domain=="exon") & (c(1:p) > intron.j)),(2*(i-1)+2)] = 1
  }
  ## Collapse
  knownDir = cbind(platDir,JDir,SVDir)
  return(knownDir)
}

#' Collect ATT/ATS directions
#'
#' @export
build_atDir = function(plat.table,inputData,Ranges,JSR.table,cutoff) {
  ## inputData = normalized data
  nexons = dim(Ranges$lRanges)[1]
  # cutoff = sqrt(qchisq(p=(1-1e-04),df=11))

  # plat.table = build_platTable(Ranges = Ranges,JSR.table = JSR.table)
  plat.baseMat = build_baseMat(plat.table)
  plat.sizeMat = sqrt(t(plat.baseMat)%*%plat.baseMat)
  normProjData = t(plat.baseMat) %*% inputData

  ## ATT
  attDir = matrix(,nrow=dim(plat.table)[1],ncol=0)
  for (i in 1:(nexons-1)) {
    intron.j = which(sapply(rownames(plat.table),FUN=function(t){strsplit(t,"[.]")[[1]][1]})==paste0("I",i))
    dir1 = dir2 = matrix(0,nrow=dim(plat.table)[1],ncol=2)
    # ATT
    if (length(intron.j)==1) {
      dir1[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) < intron.j)),1] = 1
      dir1[intron.j,2] = 1

      dir2[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) > intron.j)),1] = -1
      dir2[intron.j,2] = 1
    } else {
      dir1[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) < intron.j[1])),1] = 1
      dir1[intron.j[1],2] = 1

      dir2[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) > intron.j[1])),1] = -1
      dir2[intron.j[1],2] = 1
    }
    adj.dir1 = plat.sizeMat%*%dir1
    adj.dir2 = plat.sizeMat%*%dir2

    dir1.result = get_PO(X=t(adj.dir1)%*%normProjData,qrsc=T)
    dir2.result = get_PO(X=t(adj.dir2)%*%normProjData,qrsc=T)
    # dir1 result
    dir1.out = which(dir1.result$OS>cutoff)
    if (length(dir1.out)>0) {
      dir1.dirtmp = matrix(dir1.result$directions[,dir1.out[which((dir1.result$directions[1,dir1.out]>0) & (dir1.result$directions[2,dir1.out]>0))]],nrow=2)
      dir1.dirtmp = matrix(dir1.dirtmp[,which((dir1.dirtmp[1,]^2 > 0.01) & (dir1.dirtmp[2,]^2 > 0.01))],nrow=2)
      dir1.directions = rbind(unique(dir1.dirtmp[1,]),sqrt(1-(unique(dir1.dirtmp[1,])^2)))
      colnames(dir1.directions) = rep(paste0("ATT",i),ncol(dir1.directions))
      attDir = cbind(attDir,dir1%*%dir1.directions)
    }
    # dir2 result
    dir2.out = which(dir2.result$OS>cutoff)
    if (length(dir2.out)>0) {
      dir2.dirtmp = matrix(dir2.result$directions[,dir2.out[which((dir2.result$directions[1,dir2.out]>0) & (dir2.result$directions[2,dir2.out]>0))]],nrow=2)
      dir2.dirtmp = matrix(dir2.dirtmp[,which((dir2.dirtmp[1,]^2 > 0.01) & (dir2.dirtmp[2,]^2 > 0.01))],nrow=2)
      dir2.directions = rbind(unique(dir2.dirtmp[1,]),sqrt(1-(unique(dir2.dirtmp[1,])^2)))
      colnames(dir2.directions) = rep(paste0("ATT",i),ncol(dir2.directions))
      attDir = cbind(attDir,dir2%*%dir2.directions)
    }
    # cat(i,"\n")
  }

  ## ATS
  atsDir = matrix(,nrow=dim(plat.table)[1],ncol=0)
  for (i in 1:(nexons-1)) {
    intron.j = which(sapply(rownames(plat.table),FUN=function(t){strsplit(t,"[.]")[[1]][1]})==paste0("I",i))
    dir1 = dir2 = matrix(0,nrow=dim(plat.table)[1],ncol=2)
    # ATS
    if (length(intron.j)==1) {
      dir1[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) > intron.j)),1] = 1
      dir1[intron.j,2] = 1

      dir2[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) < intron.j)),1] = -1
      dir2[intron.j,2] = 1
    } else {
      dir1[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) > intron.j[2])),1] = 1
      dir1[intron.j[2],2] = 1

      dir2[which((plat.table$Domain=="exon") & (c(1:dim(plat.table)[1]) < intron.j[2])),1] = -1
      dir2[intron.j[2],2] = 1
    }
    adj.dir1 = plat.sizeMat%*%dir1
    adj.dir2 = plat.sizeMat%*%dir2

    dir1.result = get_PO(X=t(adj.dir1)%*%normProjData,qrsc=T)
    dir2.result = get_PO(X=t(adj.dir2)%*%normProjData,qrsc=T)
    # dir1 result
    dir1.out = which(dir1.result$OS>cutoff)
    if (length(dir1.out)>0) {
      dir1.dirtmp = matrix(dir1.result$directions[,dir1.out[which((dir1.result$directions[1,dir1.out]>0) & (dir1.result$directions[2,dir1.out]>0))]],nrow=2)
      dir1.dirtmp = matrix(dir1.dirtmp[,which((dir1.dirtmp[1,]^2 > 0.01) & (dir1.dirtmp[2,]^2 > 0.01))],nrow=2)
      dir1.directions = rbind(unique(dir1.dirtmp[1,]),sqrt(1-(unique(dir1.dirtmp[1,])^2)))
      colnames(dir1.directions) = rep(paste0("ATS",i),ncol(dir1.directions))
      atsDir = cbind(atsDir,dir1%*%dir1.directions)
    }
    # dir2 result
    dir2.out = which(dir2.result$OS>cutoff)
    if (length(dir2.out)>0) {
      dir2.dirtmp = matrix(dir2.result$directions[,dir2.out[which((dir2.result$directions[1,dir2.out]>0) & (dir2.result$directions[2,dir2.out]>0))]],nrow=2)
      dir2.dirtmp = matrix(dir2.dirtmp[,which((dir2.dirtmp[1,]^2 > 0.01) & (dir2.dirtmp[2,]^2 > 0.01))],nrow=2)
      dir2.directions = rbind(unique(dir2.dirtmp[1,]),sqrt(1-(unique(dir2.dirtmp[1,])^2)))
      colnames(dir2.directions) = rep(paste0("ATS",i),ncol(dir2.directions))
      atsDir = cbind(atsDir,dir2%*%dir2.directions)
    }
    # cat(i,"\n")
  }
  return(cbind(atsDir,attDir))
}

#' Get sparse base matrix (binary)
#'
#' @export
build_baseMat = function(plat.table) {
  plat = split_junction(as.character(plat.table$Range))
  d = max(as.numeric(plat))
  sparse.baseMat = matrix(0,ncol=dim(plat.table)[1],nrow=d)
  colnames(sparse.baseMat) = rownames(plat.table)
  for (i in 1:ncol(sparse.baseMat)) {
    tmp_bps = as.numeric(plat[,i])
    sparse.baseMat[c(tmp_bps[1]:tmp_bps[2]),i] = 1
  }
  return(sparse.baseMat)
}

#' Get plat base table for plat regions
#'
#' @export
build_platTable = function(Ranges,JSR.table) {

  plat.junctions = plat_gene(Ranges=Ranges,JSR.table=JSR.table)
  sparse.base = apply(plat.junctions,1,collapse_junction)
  sparse.base.name = rep(0,length(sparse.base))

  lRanges = Ranges$lRanges
  nexons = dim(lRanges)[1]
  for (ie in 1:nexons) {
    tmp_IDS = which((plat.junctions[,1]>=lRanges[ie,2]) & (plat.junctions[,2]<=lRanges[ie,3]))
    if (length(tmp_IDS)>1) {
      sparse.base.name[tmp_IDS] = paste(paste0("E",ie),1:length(tmp_IDS),sep=".")
    } else {
      sparse.base.name[tmp_IDS] = paste0("E",ie)
    }
    rm(tmp_IDS)
    if (ie < nexons) {
      tmp_IDS = which((plat.junctions[,1]>=lRanges[ie,3]) & (plat.junctions[,2]<=lRanges[(ie+1),2]))
      if (length(tmp_IDS)>1) {
        sparse.base.name[tmp_IDS] = paste(paste0("I",ie),1:length(tmp_IDS),sep=".")
      } else {
        sparse.base.name[tmp_IDS] = paste0("I",ie)
      }
    }
  }
  domain.name = rep("exon",length=length(sparse.base))
  domain.name[which(grepl(pattern="I",sparse.base.name))] = "intron"
  sparse.base = data.frame(Range=sparse.base,
                           Domain=domain.name,
                           Tag=sparse.base.name)
  rownames(sparse.base) = sparse.base.name
  return(sparse.base)
}

#' Get plat regions of a gene based on junctions
#'
#' @export
plat_gene = function(Ranges,JSR.table) {
  lRanges = Ranges$lRanges
  nexons = dim(lRanges)[1]
  junctions.g = matrix(as.numeric(split_junction(as.character(JSR.table$junctions.g))),ncol=2,byrow=T)
  junctions.l = matrix(as.numeric(split_junction(as.character(JSR.table$junctions.l))),ncol=2,byrow=T)

  tmp_junctions = sort(unique(c(1,c(junctions.l)[which((c(junctions.l)>0) & (c(junctions.l)<=max(lRanges)))],max(lRanges),lRanges[2:nexons,1]-1,lRanges[,3]+1)))

  junctions.pos = c(lRanges[,c(2,3)])
  for (ie in 1:nexons) {
    tmp_IDS = which((tmp_junctions>lRanges[ie,2]) & (tmp_junctions<lRanges[ie,3]))
    junctions.pos = c(junctions.pos,tmp_junctions[tmp_IDS]-1,tmp_junctions[tmp_IDS])

    if (ie < nexons) {
      if (diff(lRanges[(ie+1),c(1,2)])==0) {
        junctions.pos = c(junctions.pos,lRanges[ie,3]+1,lRanges[ie,4])
      } else {
        junctions.pos = c(junctions.pos,lRanges[ie,3]+1,lRanges[ie,4],
                          lRanges[(ie+1),1],lRanges[(ie+1),2]-1)
      }
    }
  }
  plat.junctions = matrix(sort(unique(junctions.pos)),ncol=2,byrow=T)
  return(plat.junctions)
}

#' Get skewness-adjusted PO based on A-D statistic
#'
#' @export
POrateADadj = function(x,qrsc=FALSE) {
  ADstat = ADstatWins.hy(x)
  w = max(1-(ADstat/100),0.1)
  return(w*pd.rate.hy(x,qrsc=qrsc))
}

#' Get primary PO (projection outlyingness) based on some given directions
#'
#' @export
get_POgivenB = function(X,B,qrsc=FALSE) {
  # X = d by n data matrix
  # B = d by K direction matrix
  Y = t(apply(t(B)%*%X,1,FUN= function(t) POrateADadj(t,qrsc=qrsc)))
  if (dim(Y)[1]==1) {
    return(as.vector(Y))
  } else {
    return(Y)
  }
}

#' Get sum of robust Z-scores (A-D stat adjusted) from each dimension after removing the given direction
#'
#' @export
get_resdZsum1d = function(X,direction,L1=FALSE,qrsc=TRUE) {
  resdX = X - direction%*%t(direction)%*%X
  X.stats = apply(X,1,FUN=function(t){unlist(pd.stats(t,qrsc=qrsc))})
  X.stats[3,which(grepl(pattern="I",rownames(X))==T)] = X.stats[2,which(grepl(pattern="I",rownames(X))==T)]
  resdZsum = rep(0,ncol(X))
  # sapply(resdX[1,],FUN=function(t) {compute.pd(x=t,stats=X.stats[,1])})

  ADstat = apply(X,1,ADstatWins.hy)
  ADw = sapply(ADstat,FUN=function(t){max(1-(t/100),0.1)})
  if (L1) {
    for (j in 1:nrow(resdX)) {
      resdZsum = resdZsum + ADw[j]*abs(sapply(resdX[j,],FUN=function(t) {compute.pd(x=t,stats=X.stats[,j])}))
    }
    return(resdZsum)
  } else {
    for (j in 1:nrow(resdX)) {
      resdZsum = resdZsum + (ADw[j]*sapply(resdX[j,],FUN=function(t) {compute.pd(x=t,stats=X.stats[,j])}))^2
    }
    return(sqrt(resdZsum))
  }
}

#' @export
pd.stats = function(x,qrsc=TRUE) {
  # projection depth
  if (qrsc) {
    outlier = which(abs((x-median(x))/mad(x))>qnorm(0.99))
    x0 = x[which(! 1:length(x) %in% outlier)]
    rsc=compScales(x0)
    return(list(med=rsc$med,sa=rsc$sa,sb=rsc$sb))
  } else {
    return(list(med=median(x),sa=mad(x),sb=mad(x)))
  }
}

#' @export
compute.pd = function(x,stats) {
  med = stats[1]
  sa = stats[2]
  sb = stats[3]
  y = 0
  if ((x>med) & (sa>1e-5)) {
    y = (x-med)/sa
  }
  if ((x<med) & (sb>1e-5)) {
    y = (x-med)/sb
  }
  return(y)
}

# get_resdZsum1d = function(X,direction,L1=FALSE,qrsc=FALSE) {
#   resdX = X - direction%*%t(direction)%*%X
#   resdZ = t(apply(resdX,1,FUN=function(t){POrateADadj(t,qrsc=qrsc)}))
#   if (L1) {
#     return(apply(resdZ,2,FUN=function(t){sum(abs(t))}))
#   } else {
#     return(sqrt(apply(resdZ,2,FUN=function(t){sum(t^2)})))
#   }
# }

#' Get sum of robust Z-scores (A-D stat adjusted) from each dimension after removing the given direction
#'
#' @export
get_resdZsum = function(X,B,L1=FALSE,qrsc=TRUE) {
  # X = d by n data matrix
  # B = d by n direction matrix or d-dimensional vector
  # If L1 is true (default), the L1 norm will be used to compute the sum of outlyingness. Otherwise, L2 norm will be used.
  # get_resdZsum1d = function(X,direction,qrsc=FALSE) {
  #   resdX = X - direction%*%t(direction)%*%X
  #   resdZ = t(apply(resdX,1,FUN=function(t){POrateADadj(t,qrsc=qrsc)}))
  #   if (L1) {
  #     return(apply(resdZ,2,FUN=function(t){sum(abs(t))}))
  #   } else {
  #     return(sqrt(apply(resdZ,2,FUN=function(t){sum(t^2)})))
  #   }
  # }
  if (is.null(dim(B))) {
    return(get_resdZsum1d(X=X,direction=get_unitdir(B),qrsc=qrsc))
  } else {
    unitB = apply(B,2,get_unitdir)
    return(t(apply(unitB,2,FUN=function(t){get_resdZsum1d(X=X,direction=t,qrsc=qrsc)})))
  }
}


