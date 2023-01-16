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

#' Extract the unique diplotype calls from a vector of calls
#' 
#' @param x vector of diplotype calls. i.e. c("*1|*34", *10/Amb", ...)
#' @param phased Should the output delimiter indicate that the data are phased 
#' by splitting the calls using the pipe ("|") symbol. Default TRUE. If FALSE 
#' then the "/" symbol will be used as a delimiter
#' @param ambiguous Should ambiguous calls be included in the summary. Default 
#' FALSE, removes any "Amb" calls from the final diplotype call
.uniqueCalls <- function(x, phased, ambiguous) {
    if (!is.character(x))
        x <- as.character(x)
    
    ss <- strsplit(x, "[\\|/]")
    
    h1 <- unique(unlist(lapply(ss, function(x) x[1])))
    h2 <- unique(unlist(lapply(ss, function(x) x[2])))
    
    h1 <- paste(h1, collapse = " ")
    h2 <- paste(h2, collapse = " ")
    
    if (!ambiguous)
        h1 <- gsub("Amb ", "", h1)
    h2 <- gsub("Amb ", "", h2)
    h1 <- gsub("Amb$", "", h1)
    h2 <- gsub("Amb$", "", h2)
    
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
#' and each column cooresponding to a star allele
#' @param phased logical. Default TRUE. Should the output delimiter indicate that 
#' the data was phased, i.e. should the output delimiter be "|". If FALSE the 
#' output delimter is "/"
#' @param ambiguous logical. Default FALSE. Should ambiguous ("Amb") calls be 
#' reported for each sample?
#' @return DataFrame with summarized calls for each sample
summarizeDiplotypeCalls <- function(df, phased = TRUE, ambiguous = FALSE) {
    summarized <- by(
        data = df, 
        INDICES = seq_len(nrow(df)), 
        FUN = .uniqueCalls, 
        phased = phased,
        ambiguous = ambiguous,
        simplify = FALSE
        )
    DF <- DataFrame(row.names = rownames(df), Call = as.character(summarized))
    
    return(DF)
}
