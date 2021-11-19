BuildExonsFromGTF_GeneSelection <- function(geneList = NULL, GTFfile = NULL, outputFile = NULL) {
  
  # == library setup == #
  
  req_packages <- c("ballgown", "dplyr", "GenomicRanges", "stringr", "tidyr")
  for (pack in req_packages) {
    if(pack %in% rownames(installed.packages()) == FALSE) {
      BiocManager::install(pack)
      install.packages(pack)
      suppressPackageStartupMessages(library(pack, character.only = TRUE))
    } else {
      suppressPackageStartupMessages(library(pack, character.only = TRUE))
    }
  }
  
  # == auxiliary function == #
  
  get_exons <- function(gene_id, gtfdf, gtfGenes, cols) {
    tmp.gtfdf <- gtfdf[which(gtfGenes$gene_id == gene_id), 
    ]
    tmp.bed <- tmp.gtfdf[tmp.gtfdf$feature == "exon", cols]
    exondata <- data.frame(reduce(GRanges(tmp.bed)))
    exons <- paste(exondata$seqnames[1], paste(paste(exondata$start, 
                                                     exondata$end, sep = "-"), collapse = ","), exondata$strand[1], 
                   sep = ":")
    return(exons)
  }
  
  # == subset for genes of interest + read GTF file == #
  
  if (is.null(GTFfile)) {
    stop("GTFfile must be specified.")
  }
  
  if (!is.null(geneList)) {
    # subsetting GTF for gene selection
    write.table(stringr::str_to_title(geneList), paste0(dirname(GTFfile),"/temp_genelist.txt"), 
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    system(paste0("grep -w -f ", dirname(GTFfile), "/temp_genelist.txt ", GTFfile, " > ", 
    dirname(GTFfile), "/temp.gtf"))
    
    GTFfile_temp <- paste0(dirname(GTFfile),"/temp.gtf")
  } else {
    # reading for all genes in GTF file
    GTFfile_temp <- GTFfile
  }
 
  gtfdf <- gffRead(GTFfile_temp)
  
  # get exons from only one annotation source
  ## priority: ensembl_havana > ensembl > havana
  
  suppressWarnings( # warnings due to the use of the separate function
    gtfdf <- gtfdf %>%
      tidyr::separate(attributes, c("gene_id","etc"), sep = ";", remove = FALSE) %>%
      dplyr::group_by(gene_id) %>%
      dplyr::mutate(source_priority = case_when(source == "ensembl_havana" ~ 3,
                                                source == "ensembl" ~ 2,
                                                source == "havana" ~ 1,
                                                TRUE ~ 0)) %>%
      dplyr::mutate(source_select = case_when(source_priority == max(source_priority) ~ 1,
                                              TRUE ~ 0)) %>%
      dplyr::ungroup() %>%
      dplyr::filter(source_select == 1) %>%
      dplyr::select(-gene_id, -etc, -source_priority, -source_select)
  )
  
  cols <- c("seqname", "start", "end", "strand")
  if (!"seqname" %in% colnames(gtfdf)) {
    if ("seqid" %in% colnames(gtfdf)) {
      cols <- c("seqid", "start", "end", "strand")
    }
    else {
      stop("seqid is not available")
    }
  }
  
  # == SCISSOR's method of extracting regions from GTF == #
  
  genes <- getAttributeField(gtfdf$attributes, "gene_name")
  genes <- sapply(genes, function(x) substr(x, 2, nchar(x) - 
                                              1))
  geneids <- getAttributeField(gtfdf$attributes, "gene_id")
  geneids <- sapply(geneids, function(x) substr(x, 2, nchar(x) - 
                                                  1))
  gtfGenes <- data.frame(gene_name = genes, gene_id = geneids)
  gtfGenes_distinct = distinct(gtfGenes)
  cat(paste("\t Number of genes =", dim(gtfGenes_distinct)[1]), 
      "\n")
  seqsplit <- rep(1:(round(dim(gtfGenes_distinct)[1]/1000) + 
                       1), each = 1000)
  gtfGenes_distinct <- data.frame(gtfGenes_distinct, split = seqsplit[1:dim(gtfGenes_distinct)[1]])
  gtfGenes_distinct_split <- split(x = gtfGenes_distinct, f = gtfGenes_distinct$split)
  exons_subset <- vector("list", length(gtfGenes_distinct_split))
  
  for (isp in 1:length(gtfGenes_distinct_split)) {
    gtfGenes_distinct_subset <- gtfGenes_distinct_split[[isp]]
    gtfGenes_subset <- gtfGenes[which(gtfGenes$gene_id %in% 
                                        gtfGenes_distinct_subset$gene_id), ]
    gtfdf_subset <- gtfdf[which(gtfGenes$gene_id %in% gtfGenes_distinct_subset$gene_id), 
    ]
    exons_subset[[isp]] <- data.frame(gtfGenes_distinct_subset[, 
                                                               c(1, 2)], regions = sapply(gtfGenes_distinct_subset$gene_id, 
                                                                                          function(x) get_exons(gene_id = x, gtfdf = gtfdf_subset, 
                                                                                                                gtfGenes = gtfGenes_subset, cols = cols)))
    cat(paste(c(paste(rep("=", isp), collapse = ""), paste(c(100 * 
                                                               round(isp/length(gtfGenes_distinct_split), digits = 2), 
                                                             "%")))), "\n")
  }
  all.regions <- do.call("rbind", exons_subset)
  
  # == output to file == #
  
  if (!is.null(outputFile)) {
    write.table(x = all.regions, file = outputFile, sep = "\t", 
                quote = F, row.names = F)
    cat(paste(outputFile, "has been created"), "\n")
  }
  else {
    write.table(x = all.regions, file = paste0(dirname(GTFfile), 
                                               "/BuildExonsFromGTF_GeneSelection.out"), sep = "\t", quote = F, row.names = F)
    cat(paste(paste0(dirname(GTFfile), "/BuildExonsFromGTF_GeneSelection.out"), 
              "has been created"), "\n")
  }
  
  # == remove temporary files == #

  system(paste0("rm ", dirname(GTFfile), "/temp.gtf ", dirname(GTFfile), "/temp_genelist.txt"))
  
  # == return result == #
  
  return(all.regions)
  
}
