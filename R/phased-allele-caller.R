#' Check if all haplotype ranges are in the PGx object for the given gene
#'
#' Determine callable positions by checking to see if all positions in the allele
#' definition are present in the VCF file.
#' @param x PGx object
#' @param gr GRanges object
#' @return logical vector indicating if all positions are present in the
#' sample VCF for the given haplotype.
checkCallable <- function(x, gr) {
  vcf_gr <- SummarizedExperiment::rowRanges(x)
  ov <- IRanges::subsetByOverlaps(vcf_gr, gr, type = "equal")
  if (length(ov) == length(gr)) TRUE else FALSE
}

#' Extract genotype matrix as data.table and pivot longer
#'
#' @param x PGx object
#' @return data.table with rows for every sample x position + GT
genotypeMatrixToDT <- function(x) {
  dt <- data.table::as.data.table(
    VariantAnnotation::geno(x)$GT,
    keep.rownames = "id"
    )
  dt.m <- data.table::melt(
    dt,
    id.vars = "id",
    variable.name = "sample",
    value.name = "gt",
    variable.factor = FALSE
  )
  return(dt.m)
}

#' Split the id field of a genotype data.table
#'
#' @param x data.table
splitID <- function(x) {
  x[, c("chr.vcf", "start.vcf", "end.vcf", "REF.vcf", "ALT.vcf") := tstrsplit(id, ".", fixed=TRUE)]
  return(x)
}

#' Merge the sample data.tables to the definitions
#'
#' @param x sample data.table
#' @param y haplotype definition data.table
joinSampleToDef <- function(x, y) {
  dt <- merge(x, y, by = "id", all.x = TRUE, all.y = TRUE)
  return(dt)
}

#' Clean the merged data.table
#'
#' @param x data.table after joining definition to genotype dt
#' @param gene Gene name used to determine which columns contain allele definitions
cleanJoined <- function(x, gene) {
  allele_idx <- grep(gene, names(x), value = TRUE)

  # In VCF but not in definition -- mismatch of IDs?
  x[is.na(chr.def),
    (allele_idx) := lapply(.SD, function(x) fcoalesce(x, REF.vcf)),
    by = id, .SDcols = allele_idx]
  x[is.na(chr.def), `:=`(chr.def = chr.vcf, start.def = start.vcf,
                         end.def = end.vcf, REF.def = REF.vcf,
                         ALT.def = ALT.vcf)]

  # Missing information from the definition -- fill with REF allele
  x[is.na(chr.vcf), `:=`(chr.vcf = chr.def, start.vcf = start.def,
                         end.vcf = end.def, REF.vcf = REF.def,
                         ALT.vcf = ALT.def)]

  return(x)
}

#' Convert the GT column of the cleaned data.table to nucleotide calls
#'
#' @param x data.table after cleaning
convertGTtoNucleotides <- function(x) {
  x[, c("h1", "h2") := data.table::tstrsplit(GT, "|", fixed = TRUE)]
  x[, `:=`(
    n1 = data.table::fcase(
      h1 == "0", REF.vcf,
      h1 != "0", substr(ALT.vcf, as.numeric(h1), as.numeric(h1)),
      default = NA
    ),
    n2 = data.table::fcase(
      h2 == "0", REF.vcf,
      h2 != "0", substr(ALT.vcf, as.numeric(h2), as.numeric(h2)),
      default = NA
    )
  )]

  return(x)
}

#' Hamming distance between vectors
#'
#' @param x data.table
#' @param ref reference column of the data.table
#' @param cols vector of column names containing haplotype definitions
hammingDistance <- function(x, ref, cols) {
  dt <- x[, lapply(.SD, function(x) sum(x != get(ref), na.rm = TRUE)), .SDcols = cols]
  return(dt)
}

#' Score the haplotype call based on the hamming distance value
#'
#' @param x vector of hamming distance values
#' @param mismatches logical. Return matches with a score less than or equal to this value
scoreMatches <- function(x, mismatches = 0) {
  which(x <= mismatches)
}

#' Extract the star allele string from the call string
#'
#' If the position is empty from the call string --> no exact match then assume
#' reference allele and assign the *1 allele
#'
#' @param x named vector of called star alleles
#' @param gene Gene name to remove from the star allele
#' @return character string of star allele calls for the given sample
extractStar <- function(x, gene) {
  alleles <- names(x)
  s <- gsub(gene, "", alleles)
  s <- paste(s, collapse = ",")
  if (s == "") return("*1") else return(s)
}

#' Strip suballele ids and return only unique major alleles
#'
#' @param x string of suballele calls
#' @return character string with only unique major star alleles
summarizeAlleles <- function(x) {
  s <- unlist(strsplit(x, split = ",", fixed = TRUE))
  a <- gsub("\\.[0-9]+", "", s)

  if (all(a == "*1"))
    return("*1")

  stars <- a[which(a != "*1")]
  paste(unique(stars), collapse = ",")
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
#' @param mismatches logical. Report alleles with a hamming distance (number of)
#' mismatches) less than or equal to this value. Default == 0 (perfect match).
#' @param summarize logical. Should star allele calls be summarized up to the
#' major allele level. Default TRUE. If FALSE, sub-alleles are returned for each
#' called haplotype.
#' @export
#' @return data.table with allele calls for each sample in the PGx object
callPhasedDiplotypes <- function(x, mismatches = 0, summarize = TRUE) {
  gene <- pgxGene(x)
  build <- pgxBuild(x)
  def <- switch (pgxBuild(x),
    GRCh38 = grch38_def,
    GRCh37 = grch37_def
  )
  refs <- switch (pgxBuild(x),
    GRCh38 = grch38_haplotype_grl,
    GRCh37 = grch37_haplotype_grl
  )

  # Check the sample VCF for callable haplotypes
  haplotypes <- def[Gene == gene, unique(HaplotypeName)]
  refs <- refs[haplotypes]
  callable <- sapply(refs, checkCallable, x = x, simplify = TRUE, USE.NAMES = TRUE)
  callable <- names(callable[callable])

  # Construct a definition data.table from the callable haplotypes
  dt <- def[HaplotypeName %chin% callable]
  dt.w <- dcast(dt, id ~ HaplotypeName, value.var = "VariantAllele", fill = NA)
  dt.w[, c("chr.def", "start.def", "end.def", "REF.def", "ALT.def") := tstrsplit(id, ".", fixed = TRUE)]
  dt.w[, (callable) := lapply(.SD, function(x) fcoalesce(x, REF.def)), .SDcols = callable]

  # Extract the genotype matrix from the sample VCF
  gt <- VariantAnnotation::geno(x)$GT
  gt <- data.table::as.data.table(gt, keep.rownames = "id")
  gt.m <- data.table::melt(gt, id.vars = "id", value.name = "GT", variable.name = "sample")

  # Join the callable positions onto the data
  merged <- merge(
    x = gt.m,
    y = dt.w,
    by = "id",
    all.x = TRUE,
    all.y = TRUE,
    allow.cartesian = TRUE
  )
  merged[, c("chr.vcf", "start.vcf", "end.vcf", "REF.vcf", "ALT.vcf") := tstrsplit(id, ".", fixed = TRUE)]
  merged <- unique(merged)
  cleaned <- cleanJoined(merged, gene = gene)
  converted <- convertGTtoNucleotides(cleaned)

  # Use Hamming distance as core in calling haplotypes
  hd1 <- converted[, lapply(.SD, function(x) sum(x != n1, na.rm = TRUE)), by = sample, .SDcols = callable]
  hd2 <- converted[, lapply(.SD, function(x) sum(x != n2, na.rm = TRUE)), by = sample, .SDcols = callable]
  hd1 <- as.matrix(hd1, rownames = "sample")
  hd2 <- as.matrix(hd2, rownames = "sample")
  call1 <- apply(hd1, 1, scoreMatches, mismatches = mismatches, simplify = TRUE)
  call2 <- apply(hd2, 1, scoreMatches, mismatches = mismatches, simplify = TRUE)

  # Extract star alleles from each call
  star1 <- vapply(call1, extractStar, gene = gene, FUN.VALUE = character(1))
  star2 <- vapply(call2, extractStar, gene = gene, FUN.VALUE = character(1))

  result <- data.table::data.table(sample = names(star1), H1 = star1, H2 = star2)

  if (!summarize)
    return(result)

  result[, `:=`(H1 = summarizeAlleles(H1), H2 = summarizeAlleles(H2)), by = sample]
  return(result)
}
