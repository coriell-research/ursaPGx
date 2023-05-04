#' Interface to Illumina Cyrius CYP2D6 Star Allele Caller
#'
#' @description
#' This function provides an interface to the
#' [Ilumina Cyrius CYP2D6 star allele caller](https://github.com/Illumina/Cyrius).
#' The main difference between this function and the command line version of
#' Cyrius is that logging and multi-threading are disabled in this version and 
#' this function requires the usage of a fasta reference file.
#'
#' @param files Vector of file paths to BAM/CRAM files. If BAM/CRAM files vector 
#' is named these names will be used as Sample IDs in the final output.
#' @param reference Path to the reference fasta file used in BAM/CRAM creation.
#' @param genome Reference genome. One of c("hg19", "hg37", "hg38"). Default "hg38"
#' @param output Type of output to report. One of c("simple", "verbose"). If 
#' "simple" (default) then the Cyrius "TSV" output will be returned as a 
#' DataFrame. If "verbose" then the "JSON" output will be returned as a nested 
#' list. See 'Details' for more information.
#' @param condaenv Name of the conda env containing the Cyrius dependencies. 
#' The "r-reticulate" environment is created by default when running this 
#' function for the first time. However, if you choose to use a different 
#' environment (that also contains the same Cyrius dependencies) it can be 
#' specified with this argument. Default ("r-reticulate")
#' @param ... Additional values passed to \code{reticulate::use_condaenv(condaenv, ...)}.
#' @details
#' Cyrius is a tool to genotype CYP2D6 from a whole-genome sequencing (WGS) BAM file. Cyrius uses a novel method to solve the problems caused by the high sequence similarity with the pseudogene paralog CYP2D7 and thus is able to detect all star alleles, particularly those that contain structural variants, accurately. Please refer to our [paper](https://www.nature.com/articles/s41397-020-00205-5) for details about the method.
#'
#' Cyrius has been integrated into [Illumina DRAGEN Bio-IT Platform since v3.7](https://support.illumina.com/content/dam/illumina-support/help/Illumina_DRAGEN_Bio_IT_Platform_v3_7_1000000141465/Content/SW/Informatics/Dragen/CYP2D6_Caller_fDG.htm).
#'
#' ## Interpreting the output
#'
#' ### TSV
#'
#'     | Fields in tsv     | Explanation                                                    |
#'     |:------------------|:---------------------------------------------------------------|
#'     | Sample            | Sample name                                                    |
#'     | Genotype          | Genotype call                                                  |
#'     | Filter            | Filters on the genotype call                                   |
#'
#' A genotype of "None" indicates a no-call.
#'
#'
#' There are currently four possible values for the Filter column:
#'
#'
#' - `PASS`: a passing, confident call.
#' - `More_than_one_possible_genotype`: In rare cases, Cyrius reports two possible genotypes for which it cannot distinguish one from the other. These are different sets of star alleles that result in the same set of variants that cannot be phased with short reads, e.g. *1/*46 and *43/*45. The two possible genotypes are reported together, separated by a semicolon.
#' - `Not_assigned_to_haplotypes`: In a very small portion of samples with more than two copies of CYP2D6, Cyrius calls a set of star alleles but they can be assigned to haplotypes in more than one way. Cyrius reports the star alleles joined by underscores. For example, *1_*2_*68 is reported and the actual genotype could be *1+*68/*2, *2+*68/*1 or *1+*2/*68.
#' - `LowQ_high_CN`: In rare cases, at high copy number (>=6 copies of CYP2D6), Cyrius uses less strict approximation in calling copy numbers to account for higher noise in depth and thus the genotype call could be lower confidence than usual.
#'
#' ### JSON
#'
#'     | Fields in json    | Explanation                                                    |
#'     |:------------------|:---------------------------------------------------------------|
#'     | Coverage_MAD      | Median absolute deviation of depth, measure of sample quality  |
#'     | Median_depth      | Sample median depth                                            |
#'     | Total_CN          | Total copy number of CYP2D6+CYP2D7                             |
#'     | Total_CN_raw      | Raw normalized depth of CYP2D6+CYP2D7                          |
#'     | Spacer_CN         | Copy number of CYP2D7 spacer region                            |
#'     | Spacer_CN_raw     | Raw normalized depth of CYP2D7 spacer region                   |
#'     | Variants_called   | Targeted variants called in CYP2D6                             |
#'     | CNV_group         | An identifier for the sample's CNV/fusion status               |
#'     | Variant_raw_count | Supporting reads for each variant                              |
#'     | Raw_star_allele   | Raw star allele call                                           |
#'     | d67_snp_call      | CYP2D6 copy number call at CYP2D6/7 differentiating sites      |
#'     | d67_snp_raw       | Raw CYP2D6 copy number at CYP2D6/7 differentiating sites       |
#'
#' ## Troubleshooting
#'
#' Common causes for Cyrius to produce no-calls are:
#'
#' - Low sequencing depth. We suggest a sequencing depth of 30x, which is the standard practice recommended by clinical genome sequencing.
#' - The depth of the CYP2D6/CYP2D7 region is much lower than the rest of the genome, most likely because reads are aligned to alternative contigs. If your reference genome includes alternative contigs, we suggest alt-aware alignment so that alignments to the primary assembly take precedence over alternative contigs.
#' - The majority of reads in CYP2D6/CYP2D7 region have a mapping quality of zero. This is probably due to some post-processing tools like bwa-postalt that modifies the mapQ in the BAM. We recommend using the BAM file before such post-processing steps as input to Cyrius.
#' @export
#' @md
cyrius <- function(files, reference, genome = "hg38", output = "simple", 
                   condaenv = "r-reticulate", ...) {
  stopifnot("genome must be one of c('hg19', 'hg37', hg38')" = genome %in% c("hg38", "hg37", "hg19"))
  stopifnot("output must be one of c('simple', 'verbose')" = output %in% c("simple", "verbose"))
  stopifnot("Reference fasta file does not exist" = file.exists(reference))
  
  if (!all(file.exists(files))) {
    stop("BAM/CRAM files do not exist.")
  }

  ext <- tools::file_ext(files)
  idx <- ifelse(ext == "bam", "bai", "crai")
  idx_files <- paste(files, idx, sep = ".")
  if (!all(file.exists(idx_files))) {
    stop("BAM/CRAM files are not indexed or index cannot be found in the same directory.")
  }

  # Give samples names if they do not exist
  sample_names <- names(files)
  if (is.null(sample_names)) {
    sample_names <- tools::file_path_sans_ext(basename(files))
  }

  build <- switch(genome,
    hg19 = "19",
    hg37 = "37",
    hg38 = "38"
  )

  path <- system.file("python", package = "ursaPGx")
  star_caller <- reticulate::import_from_path("star_caller", path = path, delay_load = TRUE)
  reticulate::use_condaenv(condaenv = condaenv, ...)
  
  message("Running Cyrius...")
  result <- star_caller$cyrius(
    d = reticulate::py_dict(keys = sample_names, values = files),
    genome = build,
    reference = reference,
    threads = 1L
  )
  message("Done.")

  if (output == "simple") {
    df <- S4Vectors::DataFrame(
      row.names = names(result),
      Sample = names(result),
      Genotype = vapply(result, function(x) x$Genotype, FUN.VALUE = character(1)),
      Filter = vapply(result, function(x) x$Filter, FUN.VALUE = character(1))
    )
    return(df)
  }
  return(result)
}
