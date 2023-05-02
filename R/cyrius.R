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


#' Install Cyrius Python Dependencies
#' 
#' This function will download and install the Python dependencies needed to 
#' run Cyrius from within R. 
#' @param method Installation method. By default, "auto" automatically finds a 
#' method that will work in the local environment. Change the default to force 
#' a specific installation method. Note that the "virtualenv" method is not 
#' available on Windows.
#' @param conda The path to a conda executable. Use "auto" to allow reticulate 
#' to automatically find an appropriate conda binary. See Finding Conda and 
#' \link{reticulate::conda_binary()} for more details.
#' @export
installCyrius <- function(method = "auto", conda = "auto") {
  reticulate::py_install(
      packages = c("scipy", "numpy", "pysam", "statsmodels"), 
      method = method, 
      conda = conda, 
      channel = c("conda-forge", "anaconda", "bioconda")
      )
}