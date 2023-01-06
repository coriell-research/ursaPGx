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
#' via the \code{availableHaplotypeRanges()} function.
#'
#' @param file A phased and normalized VCF file. This file must be indexed (.tbi)
#' @param gene The PGx gene to subset from VCF file. See \code{\link{availableGenes()}} for available genes.
#' @param build The genome build. One of "GRCh38" or "GRCh37".
#' @export
#' @return Object of class PGx
readPGx <- function(file, gene, build = "GRCh38") {
  stopifnot("Genome must be one of 'GRCh38' or 'GRCh37'" = build %in% c("GRCh38", "GRCh37"))
  stopifnot("Multple genes supplied to function" = length(gene) == 1)
  stopifnot("Gene not in gene list" = gene %in% availableGenes())

  build <- match.arg(build)
  ref <- switch (build,
    GRCh38 = unique(availableGeneRanges(gene, build = "GRCh38")),
    GRCh37 = unique(availableGeneRanges(gene, build = "GRCh37"))
  )
  
  # Determine the seqLevelStyle of the input file
  vcf_head <- VariantAnnotation::scanVcfHeader(file)
  vcf_style <- GenomeInfoDb::seqlevelsStyle(VariantAnnotation::reference(vcf_head))
  stopifnot("Unknown seqLevelsStyle in input VCF" = vcf_style %in% c("NCBI", "UCSC", "dbSNP", "Ensembl"))
  
  # Adjust param seqlevels if necessary
  param_style <- GenomeInfoDb::seqlevelsStyle(ref)
  if (param_style != vcf_style) {
      message("seqLevelStyles differ between reference ranges and input VCF.")
      message("Reference ranges style:", param_style)
      message("VCF style:", vcf_style)
      message("Attempting to convert", param_style, "to", vcf_style)
      GenomeInfoDb::seqlevelsStyle(ref) <- vcf_style
  }
  
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
  
  PGx(vcf_filtered, pgxGene = gene, pgxBuild = build)
}
