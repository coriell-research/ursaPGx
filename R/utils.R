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

#' Return an GRanges object containing all unique ranges for the given gene 
#' 
#' Return a GRanges object containing all unique ranges for a given PGx gene 
#' definition. This function collects all haplotype ranges for the given gene
#' and returns only ranges for unique positions across all definitions. This 
#' function is mostly intended for internal use in \code{readPGx()} where the
#' returned ranges are used to collect all positions annotated in the sample 
#' VCF.
#' 
#' @param gene Gene name
#' @param build Genome build. One of "GRCh38" or "GRCh37"
#' @return \code{GRanges} object with unique ranges for all haplotypes (star alleles)
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
  uniq <- unique(unlist(grl[[gene]]))
  
  # Return only the GRanges information -- excluding any metadata
  gr <- GenomicRanges::GRanges(
      seqnames = seqnames(uniq),
      ranges = ranges(uniq),
      strand = strand(uniq),
      seqinfo = seqinfo(uniq)
      )
  names(gr) <- NULL
  unique(gr)
}

#' Return a VRanges object of the unique ranges for the given haplotype
#' 
#' Return a VRanges object for the all ranges in the haplotype definition for 
#' the given haplotype.
#' 
#' @param haplotype Haplotype (star allele) name
#' @param build Genome build. One of "GRCh38" or "GRCh37"
#' @return \code{GRanges} object with unique ranges for the given haplotype (star allele)
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

#' Extract the unique diplotype calls from a vector of calls
#' 
#' @param x vector of diplotype calls. i.e. c("*1|*34", *10/*2", ...)
#' @param phased Should the output delimiter indicate that the data are phased 
#' by splitting the calls using the pipe ("|") symbol. Default TRUE. If FALSE 
#' then the "/" symbol will be used as a delimiter.
.uniqueCalls <- function(x, phased, ignore_nested) {
    ss <- strsplit(x, "[\\|/]")
    
    h1 <- unique(unlist(lapply(ss, function(x) x[1])))
    h2 <- unique(unlist(lapply(ss, function(x) x[2])))
    
    h1 <- paste(h1, collapse = " ")
    h2 <- paste(h2, collapse = " ")
    
    h1 <- gsub("\\*1 ", "", h1)
    h2 <- gsub("\\*1 ", "", h2)
    h1 <- gsub("\\*1$", "", h1)
    h2 <- gsub("\\*1$", "", h2)
    
    h1 <- trimws(h1, which = "both")
    h2 <- trimws(h2, which = "both")
    
    if (h1 == "") 
        h1 <- "*1"
    
    if (h2 == "")
        h2 <- "*1"
    
    delim <- "|"
    if (!phased)
        delim <- "/"
    
    summarized <- paste(h1, delim, h2)
    
    return(summarized)
}

#' Summarize a diplotype call DataFrame
#' 
#' @param df DataFrame or data.frame of diplotype calls with samples as row.names 
#' and each column corresponding to a star allele
#' @param phased logical. Default TRUE. Should the output delimiter indicate that 
#' the data was phased, i.e. should the output delimiter be "|". If FALSE the 
#' output delimter is "/"
#' @return DataFrame with summarized calls for each sample
summarizeDiplotypeCalls <- function(df, phased = TRUE) {
    by_sample <- by(
        data = df, 
        INDICES = seq_len(nrow(df)), 
        FUN = function(x) sapply(x, as.character, simplify = TRUE), 
        simplify = FALSE
        )
    calls <- sapply(by_sample, .uniqueCalls, phased = phased, simplify = TRUE)
    DF <- DataFrame(row.names = rownames(df), Call = calls)
    
    return(DF)
}
