#' Build coverage pileup based on annotations
#'
#' This function is used to get pileup data based on new transcript for the
#' purpose of the downstream outlier detection analysis.
#' @param Pileup
#' @param caseIDs
#' @param regions genomic regions formatted as chr1:1-100,200-300:+"
#' @param inputType type of intronic region contained in pileup, with choices
#'   "whole_intron", "part_intron", or "only_exon"; the first is the default.
#' @param outputType type of intronic region that will be included in output,
#'   with choices "whole_intron", "part_intron", or "only_exon"; the second is
#'   the default.
#' @keywords
#' @examples
#' regions = "chrQ:7571719-7572198,7574858-7575157,7598088-7598437:-"
#' countPileup = build_pileup(Pileup=paste0(package_path,"SCISSOR/toydata/TOY_coverage.txt"),inputType="whole_intron",caseIDs=NULL,regions=regions)
#'
#' @import BiocManager Rsamtools
#' @export
build_pileup = function(Pileup,caseIDs=NULL,regions,
                        inputType="part_intron",
                        outputType="part_intron") {
    # inputType  = "whole_intron", "part_intron", "only_exon"
    # output.ytpe = "whole_intron", "part_intron", "only_exon"

    if (missing(Pileup)) {
        stop("Pileup is missing")
    }
    if (missing(regions)) {
        stop("Regions should be specified")
    }

    strnd = strsplit(regions,":")[[1]][3]
    if (strnd=="-") {
        rawPileup = Pileup[rev(1:nrow(Pileup)),]
    } else {
        rawPileup = Pileup
    }
    if (is.null(caseIDs)) {
        caseIDs = paste0("case-",1:ncol(rawPileup),sep="")
    }
    if (inputType=="whole_intron") {

        if (outputType=="whole_intron") {
            intron.len = NULL
            ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
            covPileup = rawPileup;
        } else if (outputType=="part_intron") {
            intron.len = ceiling(len.intron.hy(regions)*0.5);
            ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
            covPileup = rawPileup[ep.new$coverage.col,] ;    #   Area to be included.
        } else if (outputType=="only_exon") {
            intron.len = 0;
            ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
            covPileup = rawPileup[ep.new$coverage.col,] ;    #   Area to be included.
        } else {
            stop(outputType," is not an option for outputType.")
        }

    } else if (inputType=="part_intron") {

        if (outputType=="whole_intron") {
            stop(outputType," is not an option when inputType=part_intron.")
        } else if (outputType=="part_intron") {
            intron.len = ceiling(len.intron.hy(regions)*0.5);
            ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
            covPileup = rawPileup ;    #   Area to be included.
            ##### If the dimension of the given pileup data is not equal to the expected one
            ##### from the intron.len calculated, should display an error. Fix this!
        } else if (outputType=="only_exon") {
            intron.len.temp = ceiling(len.intron.hy(regions)*0.5);
            ep.new.temp = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len.temp) ;

            exonic.region = c()
            for (i in 1:nrow(ep.new.temp$ep)) {
                exonic.region = c(exonic.region,(ep.new.temp$epl[i]:ep.new.temp$epr[i]));
            }
            covPileup = rawPileup[exonic.region,];
            rm(intron.len.temp,ep.new.temp,exonic.region);

            intron.len = 0;
            ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
        } else {
            stop(outputType," is not an option for outputType.")
        }

    } else if (inputType=="only_exon") {

        if (outputType=="whole_intron") {
            stop(outputType," is not an option when inputType==only_exon.")
        } else if (outputType=="part_intron") {
            stop(outputType," is not an option when inputType==only_exon.")
        } else if (outputType=="only_exon") {
            intron.len = 0;
            ep.new = find.exon.hy(regions,is.intron=TRUE,num.intron=intron.len) ;
            covPileup = rawPileup ;
        } else {
            stop(outputType," is not an option for outputType.")
        }

    } else {
        stop(inputType," is not an option for inputType.")
    }

    Ranges=get_Ranges(regions=regions,outputType=outputType)
    new.regions=Ranges$new.regions
    chr = strsplit(new.regions,":")[[1]][1]
    strtend = do.call(rbind,strsplit(strsplit(strsplit(new.regions,":")[[1]][2],",")[[1]],"-"))
    strnd = strsplit(new.regions,":")[[1]][3]
    strtend.num=matrix(as.numeric(strtend),ncol=2)
    allPos = unlist(sapply(1:nrow(strtend.num), function(x) strtend.num[x,1]:strtend.num[x,2]))

    rownames(covPileup) = allPos
    colnames(covPileup) = caseIDs
    if (strnd=="+") {
        return(covPileup)
    } else {
        return(covPileup[rev(1:nrow(covPileup)),])
    }
}
