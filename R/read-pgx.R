#' Construct row IDs for the PGx object
#'
#' The \code{readPGx()} function does not read in the ID column of the VCF by
#' default so rownames must be constructed manually. This is intentional so
#' definition ids and VCF ids can match after conversion to a genotype matrix.
#'
#' @param x PGx object
#' @return character vector of IDs
createIDs <- function(x) {
  refs <- as.character(VariantAnnotation::ref(x))
  alts <- unlist(lapply(VariantAnnotation::alt(x), function(s) paste(as.character(s), collapse = "")))

  gr <- SummarizedExperiment::rowRanges(x)
  chrs <- as.character(GenomicRanges::seqnames(gr))
  starts <- as.character(GenomicRanges::start(gr))
  ends <- as.character(GenomicRanges::end(gr))

  # Check to see that strings will not be recycled
  lens <- sapply(list(refs, alts, chrs, starts, ends), length, simplify = TRUE)
  stopifnot("All input lengths must be the same" = all(lens == lens[1]))

  # Create unique names for the rowRanges
  paste(chrs, starts, ends, refs, alts, sep = ".")
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
#' via the \code{pgxAlleleRanges()} function.
#'
#' @param file A phased and normalized VCF file. This file must be indexed (.tbi)
#' @param gene The PGx gene to subset from VCF file. See \code{\link{pgxGenes()}} for available genes.
#' @param build The genome build. One of "GRCh38" or "GRCh37".
#' @export
#' @return Object of class PGx
readPGx <- function(file, gene, build = "GRCh38") {
  stopifnot("Genome not available" = build %in% c("GRCh38", "GRCh37"))
  stopifnot("Multple genes supplied to function" = length(gene) == 1)
  stopifnot("Gene not in gene list" = gene %in% pgxGenes())

  build <- match.arg(build)
  ref <- switch (build,
    GRCh38 = unique(pgxGeneRanges(gene, build = "GRCh38")),
    GRCh37 = unique(pgxGeneRanges(gene, build = "GRCh37"))
  )
  genome <- GenomeInfoDb::genome(ref)
  tab <- Rsamtools::TabixFile(file)
  vcf <- VariantAnnotation::readVcf(tab, genome = genome, param = ref, row.names = FALSE)
  names(vcf) <- createIDs(vcf)

  PGx(vcf, pgxGene = gene, pgxBuild = build)
}
