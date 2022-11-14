#' Check if all haplotype ranges are in the PGx object for the given gene
#'
#' @param gr GRanges object
#' @param x PGx object
#' @return logical vector indicating if all positions are present in the
#' sample VCF for the given haplotype.
checkComplete <- function(gr, x) {
  ov <- IRanges::subsetByOverlaps(x, gr, type = "equal")
  if (length(ov) == length(gr))
    return(TRUE)
  return(FALSE)
}

#' Extract the genotype matrix from the PGx object and melt longer
#'
#' @import data.table
#' @param x PGx object
#' @return data.table
getGenotypeDT <- function(x) {
  gt <- VariantAnnotation::geno(x)$GT
  gt_dt <- data.table::as.data.table(gt)
  row_dt <- data.table::as.data.table(SummarizedExperiment::rowRanges(x))
  gt_dt <- cbind(row_dt, gt_dt)

  # Create unique index based on position, ref, and alt information
  gt_dt[, id := paste(seqnames, start, end, REF, paste0(unlist(ALT)), sep = ".")]

  keep_cols <- c("id", colnames(gt))
  gt_dt <- gt_dt[, .SD, .SDcols = keep_cols]
  gt_dt.m <- data.table::melt(
    gt_dt,
    id.vars = "id",
    value.name = "GT",
    variable.name = "sample",
    variable.factor = FALSE
  )
  data.table::setkey(gt_dt.m, "id")

  return(gt_dt.m)
}

#' Convert genotype calls to nucleotides
#'
#' @import data.table
#' @param x data.table
#' @return data.table
genotypeToNucleotides <- function(x) {
  dt <- copy(x)

  # Split genotype into haplotype calls
  dt[, c("h1", "h2") := tstrsplit(GT, "|", fixed = TRUE)]

  # Convert GT to nucleotides
  dt[, `:=`(n1 = data.table::fcase(
                 h1 == "0", REF,
                 h1 != "0", substr(ALT, as.numeric(h1), as.numeric(h1)),
                 default = NA),
            n2 = data.table::fcase(
                 h2 == "0", REF,
                 h2 != "0", substr(ALT, as.numeric(h2), as.numeric(h2)),
                 default = NA))]
  dt[, id := NULL]
  setcolorder(dt, c("sample", "GT", "REF", "ALT", "h1", "h2", "n1", "n2"))

  return(dt)
}

#' Score the haplotype call based on the hamming distance value
#' @param x vector of hamming distance values
#' @param exact logical. Should only 0-distance (exact matches) be returned?
scoreHaplotype <- function(x, exact = TRUE) {
  matches <- which(x == 0)
  if (!exact) {
    if (length(matches) == 0) {
      return(which.min(x))
    }
  }
  matches
}

#' Extract the star allele string from the call string
#'
#' If the position is empty from the call string --> no exact match then assume
#' reference allele and assign the *1 allele
#' @param x string of called star alleles
#' @param gene Gene name to remove from the star allele
#' @param alleles Vector of alleles that were kept in the analysis
#' @return character string of star allele calls for the given sample
extractStar <- function(x, gene, alleles) {
  s <- alleles[x]
  s <- gsub(gene, "", s)
  s <- paste(s, collapse = ",")
  if (s == "") return("*1") else return(s)
}

#' Strip suballele ids and return only unique major alleles
#' @param x string of suballele calls
#' @return character string with only unique major star alleles
summarizeAlleles <- function(x) {
  s <- unlist(strsplit(x, split = ",", fixed = TRUE))
  a <- gsub("\\.[0-9]+$", "", s)
  paste(unique(a), collapse = ",")
}

#' Call phased diplotypes (star alleles) for all samples in a PGx object
#'
#' This function attempts to call star allele diplotypes based on the phased
#' genotype data extracted from the PGx object. It first checks the sample VCF
#' for any star allele definition that is completely represented. After
#' determining 'callable' star alleles, the genotype call for each sample is
#' extracted from the sample VCF and converted to a nucleotide string for all
#' positions that are both present in the sample VCF and every callable star
#' allele. Star allele calls are then made by calculating the Hamming distance
#' between the observed string of phased haplotypes for each sample and each
#' star allele definition. By default, only star allele definitions that exactly
#' match (i.e. Hamming distance == 0) are returned. Finally, star allele
#' calls are summarized to the major allele level.
#' @import data.table
#' @param x PGx object
#' @param exact logical. Should exact matches be used to determine the star
#' allele call? default TRUE, only calls with a hamming distance == 0 are
#' returned. if FALSE then the next best match is returned. i.e. star allele
#' with the minimum hamming distance to the observed haplotype.
#' @param summarize logical. Should star allele calls be summairzed up to the
#' major allele level. Default TRUE. If FALSE, sub-alleles are returned for each
#' called haplotype.
#' @export
#' @return data.table with allele calls for each sample in the PGx object
callPhasedDiplotypes <- function(x, exact = TRUE, summarize = TRUE) {
  gene <- pgxGene(x)
  build <- pgxBuild(x)
  grl <- switch (build,
    GRCh38 = ursaPGx:::grch38_haplotype_grl,
    GRCh37 = ursaPGx:::grch37_haplotype_grl
  )
  def <- switch (build,
    GRCh38 = ursaPGx:::grch38_gene_def,
    GRCh37 = ursaPGx:::grch37_gene_def
  )
  def <- def[[gene]]

  # Get all haplotype GRanges for the desired gene
  hnames <- setdiff(names(def), c("chr", "start", "end", "REF", "ALT"))
  refs <- grl[hnames]

  # Check for callable positions
  keep <- vapply(refs, checkComplete, logical(1), x = x, USE.NAMES = TRUE)
  keep_alleles <- names(keep)[keep]

  # Keep only the callable haplotypes in the definition
  def[, id := paste(chr, start, end, REF, ALT, sep = ".")]
  keep_cols <- c("id", keep_alleles)
  def <- def[, .SD, .SDcols = keep_cols]
  data.table::setkey(def, "id")

  # Inner join the sample GTs to the definition table
  gt_dt.m <- getGenotypeDT(x)
  joined <- gt_dt.m[def, nomatch=NULL]
  joined[, c("chr", "start", "end", "REF", "ALT") := tstrsplit(id, ".", fixed=TRUE)]
  joined <- genotypeToNucleotides(joined)

  # Calculate hamming distance between observed string and definition string
  hm_n1 <- joined[, lapply(.SD, function(x) sum(x != n1)), by = sample, .SDcols = keep_alleles]
  hm_n2 <- joined[, lapply(.SD, function(x) sum(x != n2)), by = sample, .SDcols = keep_alleles]
  hm_n1 <- as.matrix(hm_n1, rownames = "sample")
  hm_n2 <- as.matrix(hm_n2, rownames = "sample")

  # Call star alleles based on hamming distance
  call1 <- apply(hm_n1, 1, scoreHaplotype, simplify = FALSE)
  call2 <- apply(hm_n2, 1, scoreHaplotype, simplify = FALSE)

  # Clean the call string by stripping the gene name
  star1 <- vapply(call1, extractStar, gene = gene, alleles = keep_alleles,
                  FUN.VALUE = character(1))
  star2 <- vapply(call2, extractStar, gene = gene, alleles = keep_alleles,
                  FUN.VALUE = character(1))

  result <- data.table::data.table(
    sample = names(star1),
    H1 = star1,
    H2 = star2
    )

  if (!summarize)
    return(result)
  result[, `:=`(H1 = summarizeAlleles(H1), H2 = summarizeAlleles(H2)), by = sample]

  return(result)
}
