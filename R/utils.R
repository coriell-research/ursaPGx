#' Return all available PGx genes
#' 
#' Return a vector of all PGx genes for which an allele definition exists.
#' 
#' @param build Genome build. One of "GRCh38" or "GRCh37"
#' @return character vector of available gene names
#' @export
#' @examples availableGenes()
availableGenes <- function(build = c("GRCh38", "GRCh37")) {
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
#' @examples availableHaplotypes()
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
#' @return \code{VRanges} object with unique ranges for all haplotypes (star alleles)
#' for the desired gene
#' @export
#' @examples availableGeneRanges("CYP2C19")
availableGeneRanges <- function(gene, build = c("GRCh38", "GRCh37")) {
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
#' @examples availableHaplotypeRanges("CYP2C19*2")
availableHaplotypeRanges <- function(haplotype, build = c("GRCh38", "GRCh37")) {
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


#' Perform diplotype calling pipeline
#' 
#' This function is a wrapper around the main pipeline steps for calling
#' diplotypes from VCF files. It wraps the following functions: (1) \code{readPGx()}
#' (2) \code{determineCallableAlleles()} (3) \code{buildReferenceDataFrame()} (4)
#' \code{convertGTtoNucleotides()} (5) \code{callPhasedDiplotypes()} and produces
#' a single DataFrame with allele calls for all samples in the VCF as output.
#' @param vcf Path to VCF file
#' @param gene The PGx gene to perform allele calling for. Must be one of \code{availableGenes()}
#' @param build The genome build. One of "GRCh38" or "GRCh37".
#' @param phased Logical value indicating whether or not the input data are 
#' phased. TRUE is the only accepted value currently.
#' @return DataFrame with sample names as rows and column as PGx gene with allele calls
#' @rdname callDiplotypes
#' @export
#' @examples
#' \dontrun{
#' # Specify the path to the VCF file
#' vcf <- "1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
#'
# Call phased diplotypes for CYP2C8
#'result <- callDiplotypes(vcf, gene = "CYP2C8", phased = TRUE)
#' 
#' }
callDiplotypes <- function(vcf, gene, build = "GRCh38", phased = TRUE) {
    stopifnot("phased = FALSE is not yet implemented" = isTRUE(phased))
    
    message("Reading in ", vcf, " as PGx object...")
    p <- readPGx(vcf, gene, build)
    message("Determining callable alleles for ", pgxGene(p), " gene...")
    p <- determineCallableAlleles(p)
    message("Building a reference from the callable alleles...")
    p <- buildReferenceDataFrame(p)
    message("Converting the genotypes to nucleotides for all samples...")
    p <- convertGTtoNucleotides(p)
    
    if (phased)
        message("Calling phased diplotypes...")
    result <- callPhasedDiplotypes(p)
    
    message("Done.")
    return(result)
}
