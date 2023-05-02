#' Interface to Illumina Cyrius CYP2D6 Star Allele Caller
#' 
#' This function provides an interface to the 
#' \href{https://github.com/Illumina/Cyrius}{Illumina Cyrius CYP2D6 Star Allele Caller}
#' @param bam_files Vector of file paths to BAM files. 
#' @param reference Path to the reference fasta file used in BAM creation.
#' @param genome Reference genome. One of c("hg19", "hg37", "hg38"). Default "hg38"
#' @param threads Number of threads to use for calling. Default (1).
#' @export
cyrius <- function(bam_files, reference, genome = "hg38", threads = 1) {
  genome <- match.arg(genome, choices = c("hg19", "hg37", "hg38"))
  message("Illumina CYP2D6 Star Caller Not Yet Implemented")
}
