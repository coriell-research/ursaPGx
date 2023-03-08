#' Show method for PGx objects
#' @rdname PGx-class
#' @export
setMethod("show", "PGx", function(object) {
  n_alleles <- length(object@pgxCallableAlleles)
  if (object@pgxCallableAlleles[1] == "") {
    n_alleles <- 0
  }

  cat(
    "Class: PGx\n",
    "PGx gene:", object@pgxGene, "\n",
    "PGx build:", object@pgxBuild, "\n",
    "Callable Alleles",
    paste0("(", n_alleles, "):"),
    head(object@pgxCallableAlleles, 3), "...\n",
    "Number of samples:", ncol(object), " \n",
    sep = " "
  )
})

#' Getter method for the pgxGene slot of a PGx object
#'
#' The 'pgxGene' slot of a PGx object contains the gene name for which the PGx
#' locations were extracted from the input VCF.
#' @export
#' @rdname pgxGene
#' @return Name of the PGx gene
setMethod("pgxGene", "PGx", function(x) x@pgxGene)

#' Getter method for the pgxBuild slot of a PGx object
#'
#' The pgxBuild slot contains the build information for which the PGx object
#' is instantiated with.
#' @export
#' @rdname pgxBuild
#' @return Name of the genome build
setMethod("pgxBuild", "PGx", function(x) x@pgxBuild)

#' Getter method for pgxCallableAlleles slot of a PGx object
#'
#' The pgxCallableAlleles slot contains a character vector of allele names for all
#' allele definitions that are completely represented in the sample VCF and are
#' thus deemed 'callable'.
#' @export
#' @rdname pgxCallableAlleles
#' @return Character vector of alleles that can be called
setMethod("pgxCallableAlleles", "PGx", function(x) x@pgxCallableAlleles)

#' Setter method for pgxCallableAlleles slot of PGx object
#'
#' @export
setReplaceMethod("pgxCallableAlleles", "PGx", function(x, value) {
  x@pgxCallableAlleles <- value
  validObject(x)
  x
})

#' Getter method for pgxReferenceDataFrame slot of a PGx object
#'
#' The pgxReferenceDataFrame contains the expanded reference strings for every
#' callable allele of the PGx object
#' @export
#' @rdname pgxReferenceDataFrame
#' @return data.frame of positions by callable haplotypes
setMethod("pgxReferenceDataFrame", "PGx", function(x) x@pgxReferenceDataFrame)

#' Setter method for pgxReferenceDataFrame slot of PGx object
#'
#' @export
setReplaceMethod("pgxReferenceDataFrame", "PGx", function(x, value) {
  x@pgxReferenceDataFrame <- value
  validObject(x)
  x
})

#' Extract a GRangesList for all defined haplotypes of the PGx object
.extractHaplotypeRanges <- function(x) {
  haplotypes <- grep(pgxGene(x), availableHaplotypes(pgxBuild(x)), value = TRUE)
  vrl <- lapply(haplotypes, availableHaplotypeRanges, build = pgxBuild(x))
  names(vrl) <- haplotypes

  return(vrl)
}

#' Determine if one set of ranges completely encompasses another
.isCallable <- function(q, s) {
  ov <- IRanges::subsetByOverlaps(q, s, type = "equal")
  if (length(unique(ov)) == length(s)) {
    return(ov)
  }
}

#' Return a vector of allele names that are able to be called for the given PGx object
#'
#' Callable alleles are alleles where all defined positions are present in the
#' sample VCF. The function will also add the vector of the callable alleles to
#' the pgxCallableAlleles slot of the PGx object
#' @export
#' @rdname determineCallableAlleles
#' @return PGx object with pgxCallableAlleles slot filled
setMethod("determineCallableAlleles", "PGx", function(x) {
  grl <- .extractHaplotypeRanges(x)
  rr <- SummarizedExperiment::rowRanges(x)
  callable <- lapply(grl, function(x) .isCallable(rr, x))
  callable <- Filter(Negate(is.null), callable)
  stopifnot("There are no callable alleles" = length(callable) >= 1)
  pgxCallableAlleles(x) <- names(callable)

  return(x)
})

#' Build the reference data frame from the callable allele positions
#'
#' The reference data.frame consists of rows for all positions present in the
#' sample VCF and columns for every callable allele haplotype. Since not all
#' positions are present in every haplotype, missing positions for each
#' haplotype are filled in with reference bases.
#' @export
#' @rdname buildReferenceDataFrame
#' @return DataFrame containing haplotype definitions
setMethod("buildReferenceDataFrame", "PGx", function(x) {
  stopifnot("There are no callable alleles" = pgxCallableAlleles(x) != "")

  grl <- lapply(pgxCallableAlleles(x), availableHaplotypeRanges, build = pgxBuild(x))
  names(grl) <- pgxCallableAlleles(x)

  rr <- SummarizedExperiment::rowRanges(x)
  df <- S4Vectors::DataFrame(REF = rr$REF, row.names = names(rr))

  for (i in seq_along(grl)) {
    idx <- SummarizedExperiment::match(grl[[i]], rr)
    alts <- VariantAnnotation::alt(grl[[i]])
    col <- names(grl)[i]
    df[idx, col] <- alts
  }
  
  coalesce <- function(x, y) ifelse(is.na(x), y, x)
  df[] <- lapply(df, function(x) coalesce(x, df$REF))
  df$REF <- NULL
  pgxReferenceDataFrame(x) <- df

  return(x)
})

#' Convert genotype calls to nucleotides
#'
#' @export
#' @rdname convertGTtoNucleotides
#' @return PGx object with geno(x)$GT converted to nucleotide representation
setMethod("convertGTtoNucleotides", "PGx", function(x) {
  REF <- VariantAnnotation::ref(x)
  ALT <- VariantAnnotation::alt(x)
  geno2 <- geno <- VariantAnnotation::geno(x)$GT

  for (i in 1:nrow(geno)) {
    geno2[i, ] <- gsub("0", as.character(REF[i]), geno[i, ])
    for (j in 1:S4Vectors::elementNROWS(ALT[i])) {
      geno2[i, ] <- gsub(
        as.character(j),
        as.character(ALT[[i]][j]),
        geno2[i, ]
      )
    }
  }
  geno(x)$GT <- geno2

  return(x)
})

#' Extract haplotype strings for a vector of genotype strings
.getHaplotypes <- function(x, sep = "|") {
  ss <- strsplit(x, sep, fixed = TRUE)
  h1 <- vapply(ss, function(x) x[1], FUN.VALUE = character(1))
  h2 <- vapply(ss, function(x) x[2], FUN.VALUE = character(1))
  l <- list(H1 = h1, H2 = h2)

  return(l)
}

#' Convert the Boolean calls to a star alleles
.boolToStar <- function(x) {
  idx <- which(x == TRUE)
  if (length(idx) > 0) {
    calls <- names(x[idx])
    calls <- regmatches(calls, regexpr("\\*[0-9]+$", calls))
    return(calls)
  }

  # Default to reference if not matches
  return("*1")
}

#' Call diplotypes from phased genotype data
#' 
#' Determine star allele calls by matching observed haplotype calls to the 
#' reference definitions for each sample in the PGx object. Diplotype calls are 
#' determined by extracting the haplotype call string for each sample and 
#' matching it against each column of the \code{pgxReferenceDataFrame(x)}. Only 
#' exact matches are reported as calls for a particular star allele. A final 
#' diplotype string is then created by concatenating haplotype calls from each
#' allele. The final output is a data.frame with rownames for each sample in the
#' PGx object and the phased diplotype calls.
#' @export
#' @rdname callPhasedDiplotypes
#' @return DataFrame of diplotype calls for each sample in the PGx object
setMethod("callPhasedDiplotypes", "PGx", function(x) {
  # Extract the converted genotype matrix from the PGx object
  GT <- geno(x)$GT

  # Split the GT matrix by Samples
  gt_by_sample <- asplit(GT, 2)

  # Extract the haplotype calls per sample
  haplotypes_by_sample <- lapply(gt_by_sample, .getHaplotypes)

  # Populate separate call matrices for each observed haplotype
  df <- pgxReferenceDataFrame(x)
  call1 <- matrix(data = FALSE, nrow = length(haplotypes_by_sample), ncol = ncol(df))
  call2 <- matrix(data = FALSE, nrow = length(haplotypes_by_sample), ncol = ncol(df))
  dimnames(call1) <- dimnames(call2) <- list(colnames(x), colnames(df))

  # Match extracted calls to definitions for each sample
  for (i in seq_along(haplotypes_by_sample)) {  # each sample
    H1 <- haplotypes_by_sample[[i]]$H1
    H2 <- haplotypes_by_sample[[i]]$H2
    for (j in seq_along(df)) {                  # each haplotype definition
      definition <- df[, j, drop = TRUE]
      if (all(H1 == definition)) {
        call1[i, j] <- TRUE
      }
      if (all(H2 == definition)) {
        call2[i, j] <- TRUE
      }
    }
  }

  # Create star allele calls for every sample
  star1 <- apply(call1, 1, .boolToStar, simplify = TRUE)
  star2 <- apply(call2, 1, .boolToStar, simplify = TRUE)

  # bind all results into a DataFrame
  result <- DataFrame(
    row.names = names(star1),
    Call = paste(star1, star2, sep = "|")
  )

  # Rename the "Call" column as the gene name
  colnames(result) <- pgxGene(x)

  return(result)
})
