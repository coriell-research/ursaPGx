#' Return all available PGx genes
#' 
#' Return a vector of all PGx genes for which an allele definition exists.
#' 
#' @param build Genome build. One of "GRCh38" or "GRCh37"
#' @return character vector of available gene names
#' @export
#' @examples pgxGenes()
availableGenes <- function(build = "GRCh38") {
  build <- match.arg(build)
  switch (build,
    GRCh38 = names(ursaPGx:::grch38_gene_grl),
    GRCh37 = names(ursaPGx:::grch37_gene_grl)
  )
}

#' Return all available PGx haplotypes (star alleles)
#'
#' Return a vector of haplotype names (star alleles) for which an allele 
#' definition exists.
#' 
#' @param build Genome build. One of "GRCh38" or "GRCh37"
#' @return character vector of available defined star alleles (haplotypes)
#' @export
#' @examples pgxHaplotypes()
availableHaplotypes <- function(build = "GRCh38") {
  build <- match.arg(build)
  switch (build,
    GRCh38 = names(ursaPGx:::grch38_haplotype_grl),
    GRCh37 = names(ursaPGx:::grch37_haplotype_grl)
  )
}

#' Return a VRanges object of the unique ranges for the given gene
#' 
#' Return a VRanges object containing all unique ranges for a given PGx gene 
#' definition. 
#' 
#' @param gene Gene name
#' @param build Genome build. One of "GRCh38" or "GRCh37"
#' @return \code{VRanges} object with unique ranges for all haplotypes (star alleles)
#' for the desired gene
#' @export
#' @examples pgxGeneRanges("CYP2C19")
availableGeneRanges <- function(gene, build = "GRCh38") {
  build <- match.arg(build)
  grl <- switch (build,
    GRCh38 = ursaPGx:::grch38_gene_grl,
    GRCh37 = ursaPGx:::grch37_gene_grl
  )
  stopifnot("Gene must be one of availableGenes()" = gene %in% availableGenes(build))
  grl[[gene]]
}

#' Return a VRanges object of the unique ranges for the given haplotype
#' 
#' Return a VRanges object for the all ranges in the haplotype definition for 
#' the given haplotype.
#' 
#' @param haplotype Haplotype (star allele) name
#' @param build Genome build. One of "GRCh38" or "GRCh37"
#' @return \code{VRanges} object with unique ranges for the given haplotype (star allele)
#' @export
#' @examples pgxHaplotypeRanges("CYP2C19_2")
availableHaplotypeRanges <- function(haplotype, build = "GRCh38") {
  build <- match.arg(build)
  grl <- switch (build,
    GRCh38 = ursaPGx:::grch38_haplotype_grl,
    GRCh37 = ursaPGx:::grch37_haplotype_grl
    )
  stopifnot("Haplotype must be one of pgxHaplotypes()" = haplotype %in% availableHaplotypes(build))
  grl[[haplotype]]
}
