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
