#' Collect ATT/ATS directions
#'
#' @export
build_atDir = function(inputData,Ranges,JSR.table,cutoff) {
  ## inputData = normalized data
  nexons = dim(Ranges$lRanges)[1]
  # cutoff = sqrt(qchisq(p=(1-1e-04),df=11))

  plat.table = build_platTable(Ranges = Ranges,JSR.table = JSR.table)
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
POrateADadj = function(x,qrsc=TRUE) {
  ADstat = ADstatWins.hy(x)
  w = max(1-(ADstat/100),0.1)
  return(w*pd.rate.hy(x,qrsc=qrsc))
}

#' Get primary PO (projection outlyingness) based on some given directions
#'
#' @export
get_POgivenB = function(X,B,qrsc=TRUE) {
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
get_resdZsum = function(X,B,L1=TRUE,qrsc=TRUE) {
  # X = d by n data matrix
  # B = d by n direction matrix or d-dimensional vector
  # If L1 is true (default), the L1 norm will be used to compute the sum of outlyingness. Otherwise, L2 norm will be used.
  get_resdZsum1d = function(X,direction,qrsc=TRUE) {
    resdX = X - direction%*%t(direction)%*%X
    resdZ = t(apply(resdX,1,FUN=function(t){POrateADadj(t,qrsc=qrsc)}))
    if (L1) {
      return(sqrt(apply(resdZ,2,FUN=function(t){sum(abs(t))})))
    } else {
      return(sqrt(apply(resdZ,2,FUN=function(t){sqrt(sum(t^2))})))
    }
  }
  if (is.null(dim(B))) {
    return(get_resdZsum1d(X=X,direction=get_unitdir(B),qrsc=qrsc))
  } else {
    unitB = apply(B,2,get_unitdir)
    return(t(apply(unitB,2,FUN=function(t){get_resdZsum1d(X=X,direction=t,qrsc=qrsc)})))
  }
}


