#' @import stringr
#'
#' @export
get_junctions = function(jsrCount,Ranges) {
  ## jsrCount = n by 1 matrix with junction-split read counts
  ## each row for each sample

  lRanges = Ranges$lRanges
  jsrCount2 = jsrCount[,1]
  caseIDs = rownames(jsrCount)
  n = length(caseIDs)

  jsrlist = sapply(jsrCount2,function(x) sapply(unlist(strsplit(x,",")),function(y) unlist(strsplit(y,":"))[1]))
  jsrsite0 = unique(unlist(jsrlist))
  jsrsite0 = jsrsite0[which(! grepl("NA",jsrsite0))] # remove NAs
  jsrsite0 = jsrsite0[which(! jsrsite0=="~")]
  jsrsite0 = jsrsite0[order(as.numeric(split_junction(jsrsite0)[1,]))]
  if (Ranges$strand=="-") {
    jsrsite0 = rev(jsrsite0)
  }
  jsrcountlist = sapply(jsrCount2, function(x) sapply(unlist(strsplit(x,",")),function(y) unlist(strsplit(y,":"))))
  extract_count = function(y,j) {
    if (length(which(y[1,]==j))>0) {
      s = as.numeric(y[2,which(y[1,]==j)])
    } else {
      s = 0
    }
  }
  jsrmat0 = sapply(jsrcountlist,function(y) sapply(jsrsite0,function(j) extract_count(y,j)))
  rownames(jsrmat0) = sapply(rownames(jsrmat0),function(x) unlist(strsplit(x,"[.]"))[1])
  jsrmat0 = jsrmat0[match(jsrsite0,rownames(jsrmat0)),]
  rownames(jsrmat0) = jsrsite0
  colnames(jsrmat0) = caseIDs
  # cat(paste("# of different splicing =",length(jsrsite0)),"\n")

  ## Junction Filtering
  ## 1. at least one sample has >= 10 reads
  mincount = 5 # minimum splice junctions required
  jsrsite = jsrsite0[which(apply(jsrmat0,1,max)>mincount)]
  jsrmat = jsrmat0[which(apply(jsrmat0,1,max)>mincount),]
  if (is.null(dim(jsrmat))) {
    jsrmat = matrix(jsrmat,nrow=1)
    colnames(jsrmat) = colnames(jsrmat0)
    rownames(jsrmat) = jsrsite
  }
  # cat(paste("# of junctions considered =",dim(jsrmat)[1]),"\n")

  ## 2. overlap with gene models
  if (dim(jsrmat)[1]==0) {
    junctions.g.name=jsrsite
    junctions.g=matrix(as.numeric(split_junction(junctions.g.name)),ncol=2,byrow=T)

    junctions.l = junctions.g
    junctions.l.name = junctions.g.name

    jsrmat.g = jsrmat
    junction.tags = as.character(NULL)
    junction.names = as.character(NULL)
  } else {
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
    if (is.null(dim(jsrmat.g))) {
      jsrmat.g = matrix(jsrmat.g,nrow=1)
      rownames(jsrmat.g) = junctions.g.name
      colnames(jsrmat.g) = colnames(jsrmat)
    }
    junction.names = paste(GeneName,"junction",1:nrow(jsrmat.g),sep="_")
    junction.tags = paste("J",1:nrow(jsrmat.g),sep="")
    rownames(jsrmat.g) = junction.names
    # jsrmat.l = jsrmat.g
    # rownames(jsrmat.l) = junctions.l.name
    # cat(paste("# of junctions considered (final) =",dim(jsrmat.g)[1]),"\n")
    # head(jsrmat.l[,1:10])
  }

  ## Annotation
  LBE.position = annotate_junction(x=junctions.l,lRanges)
  Region.tags = get_Tag(LBE.position)

  ## Junction classes
  JV.class = get_JVclass(LBE.position)
  JSR.annotation = cbind(junctions.g.name,junctions.l.name,LBE.position,JV.class,junction.tags,Region.tags)
  rownames(JSR.annotation) = junction.names
  colnames(JSR.annotation) = c("junctions.g","junctions.l","LBE.position","JV.class","JV.tag","Region.tag")

  return(list(JSR.annotation=data.frame(JSR.annotation),
              JSRmat=jsrmat.g))
}

# get_junctions = function(jsrCount,Ranges) {
#   ## jsrCount = n by 1 matrix with junction-split read counts
#   ## each row for each sample
#   require(stringr)
#   require(zoo)
#
#   lRanges = Ranges$lRanges
#   jsrCount2 = jsrCount[,1]
#   caseIDs = rownames(jsrCount)
#   n = length(caseIDs)
#   jsrsite0=c()
#   for (i in 1:n) {
#     x.split=unlist(str_split(jsrCount2[i],","))
#     jsrsite0=unique(c(jsrsite0,unlist(str_split(x.split,":"))[seq(1,(2*length(x.split)),by=2)]))
#     # cat(paste(i,"|",length(x.split)),"\n")
#   }
#   jsrsite0=sort(jsrsite0,decreasing=T)
#   # length(jsrsite0)
#
#   jsrmat0=matrix(0,ncol=n,nrow=length(jsrsite0))
#   for (i in 1:n) {
#     x=jsrCount2[i]
#     x.split=unlist(str_split(x,","))
#     x.splicing.site=unlist(str_split(x.split,":"))[seq(1,(2*length(x.split)),by=2)]
#     x.splicing.counts=unlist(str_split(x.split,":"))[seq(2,(2*length(x.split)),by=2)]
#     x.index=apply(matrix(x.splicing.site,ncol=1),1,FUN=function(x){which(jsrsite0==x)})
#     jsrmat0[x.index,i]=as.numeric(x.splicing.counts)
#   }
#   rownames(jsrmat0)=jsrsite0
#   colnames(jsrmat0)=caseIDs
#
#   # Remove NA-NA row / Remove reads from outside to outside of the gene
#   # head(jsrsite0)
#   temp = split_junction(jsrsite0)
#   if (class(temp)=="list") {
#     # remove no split-reads case.
#     if (length(which(names(temp)=="~"))>0) {
#       temp = split_junction(jsrsite0[-which(names(temp)=="~")])
#     }
#   }
#   rmj = which((temp[1,]=="NA") | (temp[2,]=="NA"))
#   rmj = c(rmj,which((apply(temp,2,max)<min(Ranges$gRanges)) | (apply(temp,2,min)>max(Ranges$gRanges))))
#   jsrsite0 = jsrsite0[which(! 1:length(jsrsite0) %in% rmj)]
#   jsrmat0 = jsrmat0[which(! 1:length(jsrsite0) %in% rmj),]
#   # cat(paste("# of different splicing =",length(jsrsite0)),"\n")
#
#   ## Junction Filtering
#   ## 1. at least one sample has >= 10 reads
#   mincount = 5 # minimum splice junctions required
#   jsrsite = jsrsite0[which(apply(jsrmat0,1,max)>mincount)]
#   jsrmat = jsrmat0[which(apply(jsrmat0,1,max)>mincount),]
#   if (is.null(dim(jsrmat))) {
#     jsrmat = matrix(jsrmat,nrow=1)
#     colnames(jsrmat) = colnames(jsrmat0)
#     rownames(jsrmat) = jsrsite
#   }
#   # cat(paste("# of junctions considered =",dim(jsrmat)[1]),"\n")
#
#   ## 2. overlap with gene models
#   if (dim(jsrmat)[1]==0) {
#     junctions.g.name=jsrsite
#     junctions.g=matrix(as.numeric(split_junction(junctions.g.name)),ncol=2,byrow=T)
#
#     junctions.l = junctions.g
#     junctions.l.name = junctions.g.name
#
#     jsrmat.g = jsrmat
#     junction.tags = as.character(NULL)
#     junction.names = as.character(NULL)
#   } else {
#     junctions.g.name=jsrsite
#     junctions.g=matrix(as.numeric(split_junction(junctions.g.name)),ncol=2,byrow=T)
#
#     junctions.l = cbind(apply(matrix(junctions.g[,1],1),1,FUN=function(x){find_region(gpos=x,Ranges=Ranges)}),apply(matrix(junctions.g[,2],1),1,FUN=function(x){find_region(gpos=x,Ranges=Ranges)}))
#     junctions.l.name = collapse_junction(junctions.l)
#
#     trow = dim(junctions.l)[1]
#     rows2exclude = which((junctions.l[,1]<1) & (junctions.l[,2]>max(lRanges)))
#
#     ## Final
#     junctions.g = junctions.g[which(! c(1:trow) %in% rows2exclude),]
#     junctions.g.name = junctions.g.name[which(! c(1:trow) %in% rows2exclude)]
#     junctions.l = junctions.l[which(! c(1:trow) %in% rows2exclude),]
#     junctions.l.name = junctions.l.name[which(! c(1:trow) %in% rows2exclude)]
#
#     jsrmat.g = jsrmat[which(! c(1:trow) %in% rows2exclude),]
#     if (is.null(dim(jsrmat.g))) {
#       jsrmat.g = matrix(jsrmat.g,nrow=1)
#       rownames(jsrmat.g) = junctions.g.name
#       colnames(jsrmat.g) = colnames(jsrmat)
#     }
#     junction.names = paste(GeneName,"junction",1:nrow(jsrmat.g),sep="_")
#     junction.tags = paste("J",1:nrow(jsrmat.g),sep="")
#     rownames(jsrmat.g) = junction.names
#     # jsrmat.l = jsrmat.g
#     # rownames(jsrmat.l) = junctions.l.name
#     # cat(paste("# of junctions considered (final) =",dim(jsrmat.g)[1]),"\n")
#     # head(jsrmat.l[,1:10])
#   }
#
#   ## Annotation
#   LBE.position = annotate_junction(x=junctions.l,lRanges)
#   Region.tags = get_Tag(LBE.position)
#
#   ## Junction classes
#   JV.class = get_JVclass(LBE.position)
#   JSR.annotation = cbind(junctions.g.name,junctions.l.name,LBE.position,JV.class,junction.tags,Region.tags)
#   rownames(JSR.annotation) = junction.names
#   colnames(JSR.annotation) = c("junctions.g","junctions.l","LBE.position","JV.class","JV.tag","Region.tag")
#
#   return(list(JSR.annotation=data.frame(JSR.annotation),
#               JSRmat=jsrmat.g))
# }

#'
#' @export
locate_region = function(tag,Ranges,JSR.table) {
  ## example:
  ## locate_region(tag="E1",Ranges=Ranges,JSR.table=JSR.table)
  ## inner function
  locate_region0 = function(tag,Ranges,JSR.table) {
    plat.table = build_platTable(Ranges,JSR.table)
    if (length(which(as.character(JSR.table$JV.tag)==tag))>0) {
      return(as.character(JSR.table$junctions.l)[which(as.character(JSR.table$JV.tag)==tag)])
    } else if ((length(which(as.character(plat.table$Tag)==tag))>0)) {
      return(as.character(plat.table$Range)[which(as.character(plat.table$Tag)==tag)])
    } else if (grepl(pattern="E",tag)) {
      num = as.numeric(strsplit(tag,"E")[[1]][2])
      return(collapse_junction(Ranges$lRanges[num,c(2,3)]))
    } else if (grepl(pattern="I",tag)) {
      num = as.numeric(strsplit(tag,"I")[[1]][2])
      return(collapse_junction(c(Ranges$lRanges[num,3]+1,Ranges$lRanges[(num+1),2]-1)))
    } else if (grepl(pattern="ATT",tag)) {
      num = as.numeric(strsplit(tag,"ATT")[[1]][2])
      return(collapse_junction(c(Ranges$lRanges[1,2],Ranges$lRanges[(num+1),2]-1)))
    } else if (grepl(pattern="ATS",tag)) {
      num = as.numeric(strsplit(tag,"ATS")[[1]][2])
      return(collapse_junction(c(Ranges$lRanges[num,3]+1,max(Ranges$lRanges))))
    } else {
      return(".")
    }
  }
  return(sapply(tag,FUN=function(t){locate_region0(tag=t,Ranges=Ranges,JSR.table=JSR.table)}))
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
    # exon.nums = as.numeric(sapply(apply(split_junction(x),1,FUN=function(z){strsplit(z,":")[[1]][1]}),
    #                               FUN=function(t){strsplit(t,"exon")[[1]][2]}))
    y = apply(split_junction(x),1,FUN=function(z){strsplit(z,":")[[1]][2]})
    if ("out" %in% y) {
      j = which(grepl(pattern="out",y)==1)
      if (! as.numeric(y[-j])==0) {
        jclass = "Cryptic_SV"
      } else {
        jclass = "Fusion_like"
      }
    } else {
      y = as.numeric(y)
      if ((y[1]<0) | (y[2]<0)) {
        jclass = "Cryptic_ES"
      } else if ((y[1]>0) | (y[2]>0)) {
        jclass = "Cryptic_IR"
      } else {
        jclass = "TBD"
      }
    }
    return(jclass)
  }
  return(sapply(LBE.position,inner_fn))
}

#'
#' @export
get_Tag = function(LBE.position) {
  # Get tags
  require(stringr)
  LBE.position = as.character(LBE.position)
  inner_fn = function(x) {
    y = apply(split_junction(x),1,FUN=function(z){strsplit(z,":")[[1]][2]})
    exon.nums = c(as.numeric(sapply(split_junction(x)[1,],FUN=function(t){strsplit(strsplit(t,":")[[1]][1],"exon")[[1]][2]})),
                  as.numeric(sapply(split_junction(x)[2,],FUN=function(t){strsplit(strsplit(t,":")[[1]][1],"exon")[[1]][2]})))
    JV.class = get_JVclass(x)
    # jtag = paste("J",i,sep="")
    jtag = "."
    if (grepl("Cryptic",JV.class)) {
      # For cryptic event
      j = which((! y== "0") & (! y== "out"))
      if (length(j)==1) {
        if (as.numeric(y[j])<0) {
          jtag = paste("E",exon.nums[j],sep="")
        } else {
          if (j==1) {
            jtag = paste("I",exon.nums[j],sep="")
          } else {
            jtag = paste("I",(exon.nums[j]-1),sep="")
          }
        }
      }
    } else if ("out" %in% y) {
      if (diff(exon.nums)==1) {
        j = which(y=="out")
        jtag = paste("E",exon.nums[j],sep="")
      }
    } else {
      # For Exon-skipping (one exon)
      if (diff(exon.nums)==2) {
        jtag = paste("E",(exon.nums[1]+1),sep="")
      }
    }
    return(jtag)
  }
  return(sapply(LBE.position,FUN=inner_fn))
}
