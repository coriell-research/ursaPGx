#' Check VCF and reference seqLevels
#'
#' Attempt to swap the reference seqLevels to match the VCF seqLevels before
#' importing the data. This step is necessary since we use ranges as a param
#' argument to the readVcf function which requires the same seqLevels in both
#' the reference ranges and the input VCF to be identical.
.checkSwapSeqLevels <- function(file, ref, build) {
  vcf_head <- VariantAnnotation::scanVcfHeader(file)
  vcf_ref <- VariantAnnotation::reference(vcf_head)
  vcf_style <- tryCatch(
    {
      style <- GenomeInfoDb::seqlevelsStyle(vcf_ref)
      if (length(style) > 1) {
        message("Multiple seqLevelStyles defined in VCF header. Using first: ", style[1])
        style <- style[1]
      }
      style
    },
    error = function(cond) {
      message("seqLevels are not defined in the VCF header.")
      style <- if (build == "GRCh37") "NCBI" else "UCSC"
      message("Attempting to use seqLevels defined by build argument: ", style)
      style
    }
  )

  stopifnot("Unknown seqLevelsStyle in input VCF" = vcf_style %in% c("NCBI", "UCSC", "dbSNP", "Ensembl"))

  # Attempt to adjust param seqlevels if necessary
  param_style <- GenomeInfoDb::seqlevelsStyle(ref)
  if (param_style != vcf_style) {
    message("seqLevelStyles differ between reference ranges and input VCF.")
    message("Reference ranges style: ", param_style)
    message("VCF style: ", vcf_style)
    message("Attempting to convert ", param_style, " to ", vcf_style)
    GenomeInfoDb::seqlevelsStyle(ref) <- vcf_style
  }

  return(ref)
}

#' Read in PGx Data from VCF File
#'
#' This function is a wrapper around \link[VariantAnnotation]{readVcf()} that
#' reads only the PGx locations into a PGx object. The user must provide an
#' indexed VCF file. The path to the indexed VCF file gets passed to
#' \link[Rsamtools]{TabixFile()} which is then used as the file argument to
#' \code{readVcf()}.
#'
#' @details The PGx locations used to subset the input VCF file are derived from
#' PharmVar VCF files available at \url{https://www.pharmvar.org/download}.
#' The "create-reference.R" file included in the data-raw directory of
#' this package describes the steps used to reproduce the construction of the
#' GRanges object used in this function. These GRanges objects are accessible
#' via the \code{availableHaplotypeRanges()} function. The \code{readPGx()} 
#' function will attempt to detect the type of seqLevels (chromosome names) used 
#' in the sample VCF (e.g. "UCSC" type - chr1, chr2, etc. or "NCBI" type, 1, 2, 
#' 3, etc.) and modify the ranges of the PharmVar reference accordingly. If the
#' VCF header does not define a reference then the function will attempt to use
#' "NCBI" type seqlevels for GRCh37 build and "UCSC" type for GRCh38 builds 
#' specified by the build argument. If multiple seqLevels are found in the VCF
#' header then the function will attempt to use the first defined level. 
#' seqLevels of the final PGx object are set to "UCSC" style for GRCh38 and 
#' "NCBI" style for GRCh37 to ensure compatibility with downstream processing.
#'
#' @param file An indexed, phased, and normalized VCF file. This file must have 
#' an index (.tbi) in the same directory.
#' @param gene The PGx gene to subset from VCF file. See 
#' \code{\link{availableGenes()}} for available genes.
#' @param build The genome build. One of "GRCh38" or "GRCh37".
#' @export
#' @return Object of class PGx
readPGx <- function(file, gene, build = c("GRCh38", "GRCh37")) {
  stopifnot("Only a single file can be used as input" = length(file) == 1)
  stopifnot("Genome must be one of 'GRCh38' or 'GRCh37'" = build %in% c("GRCh38", "GRCh37"))
  stopifnot("Multple genes supplied to function" = length(gene) == 1)
  stopifnot("Gene not in gene list" = gene %in% availableGenes())

  if (tools::file_ext(file) == "vcf") {
    stop("File extension is '.vcf'. Are you sure this is an indexed VCF file?")
  }
  
  if (gene == "CYP2D6") {
    stop("Please use `cyrius()` to perform calling for CYP2D6.")
  }

  build <- match.arg(build)
  ref <- switch(build,
    GRCh38 = availableGeneRanges(gene, build = "GRCh38"),
    GRCh37 = availableGeneRanges(gene, build = "GRCh37")
  )
  param_style <- GenomeInfoDb::seqlevelsStyle(ref)
  ref <- .checkSwapSeqLevels(file, ref, build)

  # Read in the VCF data
  genome <- GenomeInfoDb::genome(ref)
  tab <- Rsamtools::TabixFile(file)
  vcf <- VariantAnnotation::readVcf(tab, genome = genome, param = ref)

  # Perform a final check for only overlapping ranges between samples and reference
  ov <- GenomicRanges::findOverlaps(
    query = SummarizedExperiment::rowRanges(vcf),
    subject = ref,
    type = "equal"
  )
  vcf_filtered <- vcf[S4Vectors::queryHits(ov), ]

  # Convert style of sample data to reference style for downstream analysis
  GenomeInfoDb::seqlevelsStyle(vcf_filtered) <- param_style

  p <- PGx(
    vcf = vcf_filtered,
    pgxGene = gene,
    pgxBuild = build,
    pgxReferenceDataFrame = S4Vectors::DataFrame(),
    hasCallableAlleles = FALSE,
    hasReferenceDataFrame = FALSE
  )

  # Issue a warning if there are no overlapping ranges
  l <- length(SummarizedExperiment::rowRanges(p))
  if (l == 0) {
    warning("No overlapping ranges were detected between the sample VCF and the reference.")
  }

  return(p)
}
