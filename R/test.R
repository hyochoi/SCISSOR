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

  nexons = dim(Ranges$lRanges)[1]
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
                           Domain=domain.name)
  rownames(sparse.base) = sparse.base.name
  return(sparse.base)
}

#' Get plat regions of a gene based on junctions
#'
#' @export
plat_gene = function(Ranges,JSR.table) {
  lRanges = Ranges$lRanges
  nexons = dim(lRanges)[1]
  junctions.g = matrix(as.numeric(split_junction(JSR.table[,1])),ncol=2,byrow=T)
  junctions.l = matrix(as.numeric(split_junction(JSR.table[,2])),ncol=2,byrow=T)

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
