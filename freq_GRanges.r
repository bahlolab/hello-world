setwd("/wehisan/home/allstaff/z/zheng.h/my_projects/ChIP_Seq_bed_files_together//")
options(stringsAsFactors=FALSE)

# https://cran.r-project.org/web/packages/roxygen2/vignettes/rd.html

#########################################################
# formal function writing in roxygen format


#' Freqency of each GRange of a reference genome for a
#' pool of different samples
#' 
#' \code{freq_GRanges} returns the freqency of each GRange
#' in a table.
#' 
#' This is an S4 function: inputs are of class S4, output
#' is a table of integers.
#' 
#' @param ChIP_GRanges GRangesList of the reference genome
#' @param exclude_GRanges GRanges of known promoters, introns,
#' exons, etc, for which needs to be excluded
#' @return If ChIP_GRanges is a GRangesList and exclude_GRanges
#' is a GRanges, then the output will be a table of integers 
#' indicating the freqencies of each GRange in all the samples.
#' 
#' @examples
#' freq_GRanges(Ad_hmC, known_promoters)
#' freq_GRanges(enh_SCN1A, introns_exons)

freq_GRanges <- function(ChIP_GRanges, exclude_GRanges) {
  ###### check for valid inputs
  if (class(ChIP_GRanges) != "GRangesList") {
    stop("invalid GRangesList")
  }
  if (class(exclude_GRanges) != "GRanges") {
    stop("invalid exclude_GRanges")
  }
  
  ###### merge (reduce) ChIP_GRanges
  ChIP_reduced <- GRanges()
  for (i in 1:length(ChIP_GRanges)) {
    ChIP_reduced <- reduce(c(ChIP_GRanges[[i]], ChIP_reduced))
  }
  
  ###### get rid of known promoters
  if(!is.null(known_promoters)){
    enh_prom_overlap_index <- findOverlaps(ChIP_reduced, known_promoters)
    ChIP_reduced_enh_index_reverse <- unique(enh_prom_overlap_index@from)
    ChIP_reduced <- ChIP_reduced[-c(ChIP_reduced_enh_index_reverse),]
  }
  
  ###### findOverlaps(CHIP_reduced, Sample)
  Sample_histone_index <- list()
  for (i in 1:length(ChIP_GRanges)) {
    Sample_histone_index[[i]] <- (findOverlaps(ChIP_reduced, ChIP_GRanges[[i]]))@from
  }
  
  ###### count number of overlaps for each GRange in ChIP_GRanges
  tab_all_samples <- table(unlist(Sample_histone_index))
  return(tab_all_samples)
}


#########################################################
# Data processing following the function

Brain_freq <- freq_GRanges(ChIP_GRanges, known_promoters)
Brain_freq[1:10]
# 1  2  3  4  5  6  7  8  9 10 
# 3  4  4  3  5  1  1  2  7  1 
length(which(Brain1_freq>10))
# [1] 22866