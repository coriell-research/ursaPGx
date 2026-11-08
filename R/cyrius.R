#' Interface to Illumina Cyrius CYP2D6 Star Allele Caller
#'
#' @description
#' This function provides an interface to the
#' [Ilumina Cyrius CYP2D6 star allele caller](https://github.com/Illumina/Cyrius).
#' The main difference between this function and the command line version of
#' Cyrius is that logging and multi-threading are disabled in this version and
#' this function requires the usage of a fasta reference file.
#'
#' Cyrius' Python dependencies (numpy, scipy, pysam, statsmodels) are declared
#' with \code{reticulate::py_require()} when ursaPGx is loaded. On the first
#' call to \code{cyrius()} reticulate uses [uv](https://docs.astral.sh/uv/) to
#' build a cached, ephemeral Python environment satisfying these requirements
#' (an internet connection is needed the first time). No manual installation
#' is required.
#'
#' To use your own Python environment instead, select it with
#' \code{reticulate::use_condaenv()}, \code{reticulate::use_virtualenv()}, or
#' \code{reticulate::use_python()} (or set \code{RETICULATE_PYTHON}) before
#' the first call to \code{cyrius()}. That environment must provide the
#' packages listed above. Windows is not supported because pysam is not
#' available for Windows.
#'
#' @param files Vector of file paths to BAM/CRAM files. If BAM/CRAM files vector
#' is named these names will be used as Sample IDs in the final output.
#' @param reference Path to the reference fasta file used in BAM/CRAM creation.
#' Must be specified when using CRAM input. Default NULL.
#' @param genome Reference genome. One of c("hg19", "hg37", "hg38"). Default "hg38"
#' @param output Type of output to report. One of c("simple", "verbose"). If
#' "simple" (default) then the Cyrius "TSV" output will be returned as a
#' DataFrame. If "verbose" then the "JSON" output will be returned as a nested
#' list. See 'Details' for more information.
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
cyrius <- function(files, reference = NULL, genome = "hg38", output = "simple") {
  stopifnot("genome must be one of c('hg19', 'hg37', hg38')" = genome %in% c("hg38", "hg37", "hg19"))
  stopifnot("output must be one of c('simple', 'verbose')" = output %in% c("simple", "verbose"))

  if (!all(file.exists(files))) {
    stop("BAM/CRAM files do not exist.")
  }

  ext <- tolower(tools::file_ext(files))
  if (!all(ext %in% c("bam", "cram"))) {
    stop("Input files must have a .bam or .cram extension.")
  }

  # Accept the index locations htslib searches by default, e.g. for x.bam:
  # x.bam.bai, x.bai, x.bam.csi, x.csi; for x.cram: x.cram.crai, x.crai
  has_index <- vapply(seq_along(files), function(i) {
    idx_ext <- if (ext[i] == "bam") c("bai", "csi") else "crai"
    candidates <- c(
      paste(files[i], idx_ext, sep = "."),
      paste(tools::file_path_sans_ext(files[i]), idx_ext, sep = ".")
    )
    any(file.exists(candidates))
  }, FUN.VALUE = logical(1))
  if (!all(has_index)) {
    stop(
      "BAM/CRAM files are not indexed or index cannot be found in the same directory: ",
      paste(files[!has_index], collapse = ", ")
    )
  }

  # Require reference for CRAM files
  if (any(ext == "cram")) {
    if (is.null(reference)) {
      stop("Reference must be defined when using CRAM files")
    }
    stopifnot("Reference file not found" = file.exists(reference))
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

  message("Running Cyrius...")
  result <- star_caller$cyrius(
    d = reticulate::py_dict(keys = sample_names, values = files),
    genome = build,
    reference = reference,
    threads = 1L
  )
  message("Done.")

  if (output == "simple") {
    tryCatch(
      {
        df <- S4Vectors::DataFrame(
          row.names = names(result),
          Sample = names(result),
          Genotype = vapply(result, function(x) x$Genotype, FUN.VALUE = character(1)),
          Filter = vapply(result, function(x) x$Filter, FUN.VALUE = character(1))
        )
        return(df)
      },
      error = function(cond) {
        message("An error occurred when extracting json results for at least one of the samples:")
        message(conditionMessage(cond))
        message("")
        message("Returning verbose output instead of simple output.")
        return(result)
      }
    )
  }
  return(result)
}
