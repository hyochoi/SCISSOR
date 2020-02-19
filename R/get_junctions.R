#'
#' @export
get_junctions = function(jsrCount,Ranges) {
  ## jsrCount = n by 1 matrix with junction-split read counts
  ## each row for each sample
  require(stringr)
  require(zoo)

  lRanges = Ranges$lRanges
  jsrCount2 = jsrCount[,1]
  caseIDs = rownames(jsrCount)
  n = length(jsrCount2)
  jsrsite0=c()
  for (i in 1:n) {
    x=jsrCount2[i]
    x.split=unlist(str_split(x,","))
    jsrsite0=unique(c(jsrsite0,unlist(str_split(x.split,":"))[seq(1,(2*length(x.split)),by=2)]))
    # cat(paste(i,"|",length(x.split)),"\n")
  }
  jsrsite0=sort(jsrsite0,decreasing=T)
  # length(jsrsite0)

  jsrmat0=matrix(0,ncol=n,nrow=length(jsrsite0))
  for (i in 1:n) {
    x=jsrCount2[i]
    x.split=unlist(str_split(x,","))
    x.splicing.site=unlist(str_split(x.split,":"))[seq(1,(2*length(x.split)),by=2)]
    x.splicing.counts=unlist(str_split(x.split,":"))[seq(2,(2*length(x.split)),by=2)]
    x.index=apply(matrix(x.splicing.site,ncol=1),1,FUN=function(x){which(jsrsite0==x)})
    jsrmat0[x.index,i]=as.numeric(x.splicing.counts)
  }
  rownames(jsrmat0)=jsrsite0
  colnames(jsrmat0)=caseIDs

  # Remove NA-NA row
  # head(jsrsite0)
  temp = split_junction(jsrsite0)
  rmj = which((temp[1,]=="NA") | (temp[2,]=="NA"))
  jsrsite0 = jsrsite0[-rmj]
  jsrmat0 = jsrmat0[-rmj,]
  # cat(paste("# of different splicing =",length(jsrsite0)),"\n")

  ## Junction Filtering
  ## 1. at least one sample has >= 10 reads
  mincount = 5 # minimum splice junctions required
  jsrsite = jsrsite0[which(apply(jsrmat0,1,max)>mincount)]
  jsrmat = jsrmat0[which(apply(jsrmat0,1,max)>mincount),]
  # cat(paste("# of junctions considered =",dim(jsrmat)[1]),"\n")

  ## 2. overlap with gene models
  junctions.g.name=jsrsite
  junctions.g=matrix(as.numeric(split_junction(junctions.g.name)),ncol=2,byrow=T)

  junctions.l = cbind(apply(matrix(junctions.g[,1],1),1,FUN=function(x){find_region(gpos=x,Ranges=Ranges)}),apply(matrix(junctions.g[,2],1),1,FUN=function(x){find_region(gpos=x,Ranges=Ranges)}))
  junctions.l.name = collapse_junction(junctions.l)

  trow = dim(junctions.l)[1]
  rows2exclude = which((junctions.l[,1]<1) & (junctions.l[,2]>max(lRanges)))

  ## Final
  junctions.g = junctions.g[which(! c(1:trow) %in% rows2exclude),]
  junctions.g.name = junctions.g.name[which(! c(1:trow) %in% rows2exclude)]
  junctions.l = junctions.l[which(! c(1:trow) %in% rows2exclude),]
  junctions.l.name = junctions.l.name[which(! c(1:trow) %in% rows2exclude)]

  jsrmat.g = jsrmat[which(! c(1:trow) %in% rows2exclude),]
  junction.names = paste(GeneName,"junction",1:nrow(jsrmat.g),sep="_")
  rownames(jsrmat.g) = junction.names
  # jsrmat.l = jsrmat.g
  # rownames(jsrmat.l) = junctions.l.name
  # cat(paste("# of junctions considered (final) =",dim(jsrmat.g)[1]),"\n")
  # head(jsrmat.l[,1:10])

  ## Annotation
  LBE.position = annotate_junction(x=junctions.l,lRanges)
  cbind(junctions.g.name,junctions.l.name,LBE.position)

  ## Junction classes
  JV.class = get_JVclass(LBE.position)
  JSR.annotation = cbind(junctions.g.name,junctions.l.name,LBE.position,JVclass)
  rownames(JSR.annotation) = junction.names
  colnames(JSR.annotation) = c("junctions.g","junctions.l","LBE.position","JV.class")
  return(list(JSR.annotation=JSR.annotation,
              JSRmat=jsrmat.g))
}

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
  nexons = nrow(lRanges)
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
get_JVclass = function(LBE.position) {
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
