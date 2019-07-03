#' Find intron length to be inserted
#'
#' @export
len.intron.hy = function(exon){
  # 2nd version. (undated on 4/2/2017)
  # Input: exon info.
  # Produces the length of introns whose total length is equal to the total length of exons
  a <- exon ;
  b <- strsplit(a,":")[[1]][2] ;
  c <- gsub("-",",",b) ;
  d <- as.numeric(strsplit(c,",")[[1]]) ;
  epl <- d[seq(1,length(d),by=2)] ;
  epr <- d[seq(2,length(d),by=2)] ;
  ep <- cbind(epl,epr) ;

  n.ex = nrow(ep);
  if (n.ex==1){
    return(0)
  } else {
    n.int = nrow(ep)-1;
    len.int = epl[2:n.ex]-epr[1:(n.ex-1)] - 1 ; # Intron lengths between exons
    srt.len.int = sort(len.int) ; # sorted intron lengths between exons
    ttl.ex = sum(epr - epl + 1) ; # total length of exons
    mean.int = floor(ttl.ex/n.int) ; # Expected mean of intron lengths

    if (srt.len.int[1]>mean.int){
      output = mean.int ;
    } else if (srt.len.int[n.int] <= mean.int){
      output = srt.len.int[n.int] ;
    } else if (sum(len.int)<ttl.ex){
      output = srt.len.int[n.int];
    } else {
      jp = 0; jc = 1;
      while (jc > jp){
        jp = jc ;
        crr.ttl.int = ttl.ex - sum(srt.len.int[1:jc]) ; # total lengths of introns at current step
        crr.mean.int = floor(crr.ttl.int/(n.int-jc)) ;

        if (srt.len.int[jc+1]>crr.mean.int){
          output = crr.mean.int ;
          jc = jp ;
        } else {
          jc = jp + 1 ;
        }
      }
    }
    return(output) ;
  }
}
