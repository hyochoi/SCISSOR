#'
#' @export
split_junction = function(x) {
  require(stringr)
  inner_fun = function(x) {strsplit(x,"~")[[1]]}
  return(sapply(x,inner_fun))
}

#'
#' @export
collapse_junction = function(x) {
  # each row for junction positions
  require(stringr)
  inner_fun = function(x) {paste(x,collapse="~")}
  if (is.null(dim(x))) {
    return(inner_fun(x))
  } else {
    return(apply(x,1,inner_fun))
  }
}

#'
#' @export
annotate_junction = function(x,lRanges) {
  # Each row for each junction (in locus)
  require(stringr)
  inner_fn = function(x,lRanges) {
    ie = which((lRanges[,4]-x[1])*(lRanges[,1]-x[1])<=0)
    if (length(ie)==0) {
      y1 = "exon1:out"
    } else if (x[1]>=lRanges[ie,3]) {
      lbe.diff = x[1] - lRanges[ie,3]
      y1 = paste0("exon",ie,":",lbe.diff)
    } else if (x[1]>=lRanges[ie,2]) {
      lbe.diff = x[1] - lRanges[ie,3]
      y1 = paste0("exon",ie,":",lbe.diff)
    } else if (0<x[1]) {
      lbe.diff = x[1] - lRanges[(ie-1),3]
      y1 = paste0("exon",(ie-1),":",lbe.diff)
    } else {
      y1 = "exon1:out"
    }
    rm(ie)
    ie = which((lRanges[,4]-x[2])*(lRanges[,1]-x[2])<=0)
    if (length(ie)==0) {
      y2 = paste0("exon",nexons,":out")
    } else if (x[2]<=lRanges[ie,2]) {
      lbe.diff = lRanges[ie,2] - x[2]
      y2 = paste0("exon",ie,":",lbe.diff)
    } else if (x[2]<=lRanges[ie,3]) {
      lbe.diff = lRanges[ie,2] - x[2]
      y2 = paste0("exon",ie,":",lbe.diff)
    } else if (x[2]<=max(lRanges)) {
      lbe.diff = lRanges[(ie+1),2] - x[2]
      y2 = paste0("exon",(ie+1),":",lbe.diff)
    } else {
      y2 = paste0("exon",nexons,":out")
    }
    return(collapse_junction(c(y1,y2)))
  }
  if (is.null(dim(x))) {
    return(inner_fn(x=x,lRanges=lRanges))
  } else {
    return(apply(x,1,FUN=function(y){inner_fn(x=y,lRanges=lRanges)}))
  }
}

#'
#' @export
JV_class = function(LBE.position) {
  # Get junction variants classes
  require(stringr)
  inner_fn = function(x) {
    y = apply(split_junction(x),1,FUN=function(z){strsplit(z,":")[[1]][2]})
    if ("out" %in% y) {
      jclass = "SV"
    } else {
      y = as.numeric(y)
      if ((y[1]<0) | (y[2]<0)) {
        jclass = "Cryptic_ES"
      } else if ((y[1]>0) | (y[2]>0)) {
        jclass = "Cryptic_IR"
      } else {
        jclass = "None"
      }
    }
    return(jclass)
  }
  return(sapply(LBE.position,inner_fn))
}
