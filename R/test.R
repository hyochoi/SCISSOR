#' Get flattened regions of a gene based on junctions
#'
#' @export
flatten_gene = function(Ranges,JSR.table) {
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
  flat.junctions.l = sort(unique(junctions.pos))
  return(flat.junctions.l)
}

#' Get primary PO (projection outlyingness) based on some given directions
#'
#' @export
get_POgivenB = function(X,B,qrsc=TRUE) {
  # X = d by n data matrix
  # B = d by K direction matrix
  Y = t(apply(t(B)%*%X,1,FUN= function(t) pd.rate.hy(t,qrsc=qrsc)))
  if (dim(Y)[1]==1) {
    return(as.vector(Y))
  } else {
    return(Y)
  }
}

#' Get a weighted sum of robust Z-scores from each dimension after removing the given direction
#'
#' @export
get_resdZsum = function(X,B,weight=NULL,qrsc=TRUE) {
  # X = d by n data matrix
  # B = d by n direction matrix or d-dimensional vector
  # weight = d-dimensional vector, weights for each dimension, used to calculate weighted sum
  if (is.null(weight)) {
    weight=rep(1,dim(X)[1])
  }
  get_resdZsum1d = function(X,direction,qrsc=TRUE) {
    resdX = X - direction%*%t(direction)%*%X
    resdZ = t(apply(resdX,1,FUN=function(t){pd.rate.hy(t,qrsc=qrsc)}))
    return(sqrt(apply(resdZ,2,FUN=function(t){sqrt(sum(weight*(t^2)))})))
  }
  if (is.null(dim(B))) {
    return(get_resdZsum1d(X=X,direction=get_unitdir(B),qrsc=qrsc))
  } else {
    unitB = apply(B,2,get_unitdir)
    return(t(apply(unitB,2,FUN=function(t){get_resdZsum1d(X=X,direction=t,qrsc=qrsc)})))
  }
}
