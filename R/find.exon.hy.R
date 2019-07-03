#' Finding exons and introns
#'
#' This function is used to find exons and introns that will be used
#' in the downstream analysis
#' @export
find.exon.hy <- function(exon,is.intron=FALSE,num.intron=NULL){
  ## % First version: find.exon.hy.2; Second version: find.exon.hy.3;
  ## % Last updated: 1/14/2016
  ## % Find the locatons of the exon boundaries in plot
  ## % If there are only exons, then is.intron=FALSE
  ## % If the coverage file contains some intron, is.intron=TRUE and the number of intron is indicated as intron.num
  ## % If intron.num=NULL, whole introns are contained.
  a <- exon ;
  b <- strsplit(a,":")[[1]][2] ;
  c <- gsub("-",",",b) ;
  d <- as.numeric(strsplit(c,",")[[1]]) ;
  epl <- d[seq(1,length(d),by=2)] ;
  epr <- d[seq(2,length(d),by=2)] ;
  ep <- cbind(epl,epr) ;
  if (!is.intron){
    ep2 <- ep - ep[1,1] + 1;
    ep3 <- ep - ep[,1] + 1 ;
    if (nrow(ep3)>1){ new.ep <- ep3 + c(0,cumsum(ep3[1:(nrow(ep3)-1),2])) } else { new.ep <- ep3 }
    coverage.col <- c();
    for (iep in 1:nrow(ep2)){
      coverage.col <- c(coverage.col, ep2[iep,1]:ep2[iep,2])
    }
    return(list(coverage.col=coverage.col, epl=new.ep[,1], epr=new.ep[,2]))
  } else {
    ep2 <- ep - ep[1,1] + 1;
    if (is.null(num.intron)){
      nep2 = nrow(ep2)
      if (nep2>1) {
        ep2 <- cbind(ep2[,1],ep2,c(ep2[2:nep2,1]-1,ep2[nep2,2]))
      } else {
        ep2 <- matrix(c(1,ep2,ep2[2]),ncol=4)
      }
      new.ep <- ep2;
      # return(list(epl=new.ep[,1],epr=new.ep[,2]))
    } else {
      ep2 <- cbind(ep2[,1]-num.intron,ep2,ep2[,2]+num.intron); ep2[1,1] <- 1;
      if (nrow(ep2) > 1){
        overlap.id <- which(ep2[1:(nrow(ep2)-1),4]>ep2[2:nrow(ep2),1]) ;
        if (length(overlap.id)>0){
          ep2[overlap.id,4] <- ep2[overlap.id,3];
          ep2[(overlap.id+1),1] <- ep2[overlap.id,3] + 1;
        }
        ep2[nrow(ep2),4] <- ep2[nrow(ep2),3];
        ep3 <- ep2 - ep2[,1] + 1;
        new.ep <- ep3 + c(0,cumsum(ep3[1:(nrow(ep3)-1),4]));
        overlap.id.2 <- which(new.ep[1:(nrow(new.ep)-1),4]>new.ep[2:nrow(new.ep),2]) ;
        if (length(overlap.id.2)>0){
          print(paste("!!!! Warning: Some exons and introns are overlapped!!!!  Overlapped at:  ",overlap.id.2))
        }
      } else {
        ep2[nrow(ep2),4] <- ep2[nrow(ep2),3];
        new.ep <- ep2
      }
    }
    #	Positions will be used for introns and exons
    coverage.col <- c();
    for (iep in 1:nrow(ep2)){
      coverage.col <- c(coverage.col, ep2[iep,1]:ep2[iep,4])
    }
    return(list(coverage.col=coverage.col,epl=new.ep[,2],epr=new.ep[,3],ep=new.ep));
    #	coverage.col = bp positions that will be used among the whole transcript (exon+intron)
    #	epl, epr = left and right exon positions in coverage.col

  }
}
