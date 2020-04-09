#' @export
get_SRcount = function(BAMfiles,caseIDs=NULL,regions=NULL,print.proc=FALSE) {
  if (missing(BAMfiles)) {
    stop("BAM file path is missing")
  }
  if (is.null(caseIDs)) {
    caseIDs = 1:length(BAMfiles)
  } else {
    if (! length(BAMfiles)==length(caseIDs)) {
      stop("BAM files and case IDs have different lengths")
    }
  }

  Ranges = get_Ranges(regions=regions)
  gene.locus = paste(Ranges$chr,paste((min(Ranges$gRanges)-100),(max(Ranges$gRanges)+100),sep="-"),Ranges$strand,sep=":")

  splicedReads.count=rep(0,length(caseIDs))
  if (print.proc) {
    for (sj in 1:length(caseIDs)) {
      bam.jsr = extract_splicedReads(BAM=BAMfiles[sj],gene.locus=gene.locus)
      splicedReads.count[sj] = count_splicedReads(splicedReads=bam.jsr,regions=regions)
      cat(paste0(sj,"|",caseIDs[sj]),"\n")
      rm(bam.jsr)
    }
  } else {
    for (sj in 1:length(caseIDs)) {
      bam.jsr = extract_splicedReads(BAM=BAMfiles[sj],gene.locus=gene.locus)
      splicedReads.count[sj] = count_splicedReads(splicedReads=bam.jsr,regions=regions)
      rm(bam.jsr)
    }
  }
  splicedReads.count = data.frame(caseIDs=caseIDs,junction.count=splicedReads.count)
  return(splicedReads.count)
}

## example
# manifest = read.table(paste0(base,"HNSC_RNA_manifest_datastore.txt"),
#                       header=F,colClasses="character")
# bamfiles = manifest[,2]
# BAM = bamfiles[1]
# gene.locus = "chr4:187508000-187650000:-"
# bam.jsr = extract_splicedReads(BAM=BAM,gene.locus=gene.locus)
# p0 = c("qname","flag","rname","pos","mapq","cigar","mrnm","mpos","isize","seq")
# p1 = c("qname","flag","rname","pos","mapq","cigar")
#' @export
extract_splicedReads = function(BAM,gene.locus,p1=c("qname","flag","rname","pos","mapq","cigar")) {
  require(Rsamtools)

  bf = BamFile(BAM)
  chr = strsplit(gene.locus,":")[[1]][1]
  strtend = do.call(rbind,strsplit(strsplit(strsplit(gene.locus,":")[[1]][2],",")[[1]],"-"))
  strnd = strsplit(gene.locus,":")[[1]][3]
  df = GRanges(chr,IRanges(start=as.numeric(strtend[,1]), end=as.numeric(strtend[,2])),strnd)
  s_param = ScanBamParam(which=df, what=p1)

  res = scanBam(bf, param=s_param)[[1]]
  output = data.frame(matrix(unlist(res),ncol=length(res)))
  colnames(output) = p1
  output = output[which(grepl(pattern="N",output$cigar)==T),]
  return(output[which(as.numeric(as.character(output$mapq))>30),])
}

#' @export
count_splicedReads=function(splicedReads,regions) {
  require("stringr")
  seqdir=unlist(str_split(regions,":"))[3]
  splicedReads = splicedReads[which(! as.character(splicedReads$cigar)==""),]

  inner_fn=function(rinfo) {
    rpos = as.numeric(as.character(rinfo[1]))
    rcigar = unlist(str_split(as.character(rinfo[2]),"[M,N]"))
    if (length(rcigar)==4) {
      rcigar = as.numeric(rcigar)
      y=c((rpos+rcigar[1]-1),
          rcigar[1:3],
          (rpos+rcigar[1]+rcigar[2]))
      # y = c(j1.pos,j1.len,split.len,j2.len,j2.pos)
    } else if (length(rcigar)==6) {
      rcigar = as.numeric(rcigar)
      y1=c((rpos+rcigar[1]-1),
           rcigar[1:3],
           (rpos+rcigar[1]+rcigar[2]))
      y2=c((rpos+sum(rcigar[1:3])-1),
           rcigar[3:5],
           (rpos+sum(rcigar[1:4])))
      y=t(rbind(y1,y2))
    } else {
      y=rep(NA,5)
    }
    return(y)
  }
  splicedReads2 = cbind(as.character(splicedReads$pos),as.character(splicedReads$cigar))
  # temp_out=apply(subset(splicedReads, select=c(pos,cigar)),1,inner_fn)
  temp_out=apply(splicedReads2,1,inner_fn)
  splice.info = matrix(unlist(temp_out),ncol=5,byrow=T)

  # Get unique splicing locations
  if (length(which(splice.info[,1]>splice.info[,5]))>0) {
    stop("bam slice should be resorted...")
  }

  if (seqdir=="+") {
    splicing.set=sort(unique(paste0(splice.info[,1],"~",splice.info[,5])),decreasing=F)
  } else {
    splicing.set=sort(unique(paste0(splice.info[,5],"~",splice.info[,1])),decreasing=T)
  }

  # Count spliced reads
  splice.count=rep(0,length(splicing.set))
  for (i in 1:length(splicing.set)) {
    splice_site=unlist(str_split(splicing.set[i],"~"))
    temp.count=length(which((splice.info[,1]==min(splice_site)) & (splice.info[,5]==max(splice_site))))
    splice.count[i]=paste0(splicing.set[i],":",temp.count)
  }
  return(paste(splice.count,collapse=","))
}
