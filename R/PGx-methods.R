#' Show method for PGx objects
#' @rdname PGx-class
#' @export
setMethod("show", "PGx", function(object) {
  cat(
    "Class: PGx\n",
    "PGx gene:", object@pgxGene, "\n",
    "PGx build:", object@pgxBuild, "\n",
    "Callable Alleles",
    paste0("(", length(object@callableAlleles), "):"),
    head(object@callableAlleles, 3), "...\n",
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

#' Getter method for callableAlleles slot of a PGx object
#'
#' The callableAlleles slot contains a character vector of allele names for all
#' allele definitions that are completely represented in the sample VCF and are
#' thus deemed 'callable'.
#' @export
#' @rdname callableAlleles
#' @return Character vector of alleles that can be called
setMethod("callableAlleles", "PGx", function(x) x@callableAlleles)

#' Setter method for callableAlleles slot of PGx object
#'
#' @export
setReplaceMethod("callableAlleles", "PGx", function(x, value) {
  x@callableAlleles <- value
  validObject(x)
  x
})

#' Extract a GRangesList for all defined haplotypes of the PGx object
#'
#' @export
#' @rdname extractHaplotypeRanges
#' @return GRangesList of haplotype definitions
setMethod("extractHaplotypeRanges", "PGx", function(x) {
  haplotypes <- grep(pgxGene(x), availableHaplotypes(pgxBuild(x)), value = TRUE)
  vrl <- lapply(haplotypes, availableHaplotypeRanges, build = pgxBuild(x))
  names(vrl) <- haplotypes

  return(vrl)
})

#' Determine if one set of ranges completely encompasses another
.isCallable <- function(q, s) {
  ov <- IRanges::subsetByOverlaps(q, s, type = "equal")
  if (length(ov) == length(s)) {
    return(ov)
  }
}

#' Return a vector of allele names that able to be called for the given PGx object
#'
#' Callable alleles are alleles where all defined positions are present in the
#' sample VCF. The function will also add the vector of the callable alleles to
#' the callableAlleles slot of the PGx object
#' @export
#' @rdname getCallableAlleles
#' @return PGx object with callableAlleles slot filled
setMethod("getCallableAlleles", "PGx", function(x) {
  grl <- extractHaplotypeRanges(x)
  q <- SummarizedExperiment::rowRanges(x)
  callable <- lapply(grl, function(x) .isCallable(q, s = x))
  callable <- Filter(Negate(is.null), callable)
  stopifnot("There are no callable alleles" = length(callable) >= 1)
  callableAlleles(x) <- names(callable)

  return(x)
})

#' Taken from:
#' https://github.com/Bioconductor/VariantAnnotation/blob/54d2de6c7a1aa76f404076db79077cd05114edf2/R/methods-readVcf.R
.geno2geno <- function(lst, ALT = NULL, REF = NULL, GT = NULL) {
  if (is.null(ALT) && is.null(REF) && is.null(GT)) {
    ALT <- lst$ALT
    REF <- as.character(lst$REF, use.names = FALSE)
    GT <- lst$GENO$GT
  }
  res <- GT
  ## ignore records with GT ".|."
  if (any(missing <- grepl(".", GT, fixed = TRUE))) {
    GT[missing] <- ".|."
  }
  phasing <- rep("|", length(GT))
  phasing[grepl("/", GT, fixed = TRUE)] <- "/"

  ## replace
  GTstr <- strsplit(as.vector(GT), "[|,/]")
  if (any(elementNROWS(GTstr) != 2)) {
    stop("only diploid variants are supported")
  }
  GTmat <- matrix(unlist(GTstr), ncol = 2, byrow = TRUE)
  GTA <- suppressWarnings(as.numeric(GTmat[, 1]))
  GTB <- suppressWarnings(as.numeric(GTmat[, 2]))

  REFcs <- cumsum(elementNROWS(REF))
  ALTcs <- cumsum(elementNROWS(ALT))
  cs <- REFcs + c(0, head(ALTcs, -1))
  offset <- rep(cs, ncol(res))
  alleles <- unlist(rbind(REF, ALT), use.names = FALSE)

  alleleA <- alleles[offset + GTA]
  alleleB <- alleles[offset + GTB]
  if (any(missing)) {
    res[!missing] <- paste0(
      alleleA[!missing],
      phasing[!missing],
      alleleB[!missing]
    )
  } else {
    res[] <- paste0(alleleA, phasing, alleleB)
  }
  res
}

#' Converts the 'GT' genotype code in a PGx object to nucleotides
#'
#' Modified from: #' https://github.com/Bioconductor/VariantAnnotation/blob/54d2de6c7a1aa76f404076db79077cd05114edf2/R/methods-VCF-class.R
#' This method differs in that it uses the REF and ALT columns from the PGx
#' haplotype definitions instead of the ALT and REF from the sample VCF and
#' requires that a specific callable allele be used. This method returns a
#' PGx object that is subsetted by the defined allele, i.e. the input PGx object
#' is not the same dimensions as the returned PGx object
#' @export
#' @rdname pgxGenotypeCodesToNucleotides
#' @return subsetted PGx object with geno(PGx)$GT converted to nucleotides
setMethod("pgxGenotypeCodesToNucleotides", "PGx", function(x, allele, ...) {
  # Extract the PGx definition for the given allele
  if (!allele %in% callableAlleles(x)) {
    stop("allele must be one of callableAlleles(x)")
  }
  def_gr <- availableHaplotypeRanges(allele, build = pgxBuild(x))

  # Subset the original PGx object for only the ranges in the definition
  p <- IRanges::subsetByOverlaps(x, def_gr, type = "equal")
  
  # Reorder the definition ranges to match the PGx ranges
  def_gr <- def_gr[GenomicRanges::match(SummarizedExperiment::rowRanges(p), def_gr), ]
  def_alt <- Biostrings::DNAStringSetList(as.list(def_gr$`Variant Allele`))
  def_ref <- def_gr$`Reference Allele`

  # Convert
  GT <- VariantAnnotation::geno(p)$GT
  ALT <- as.list(S4Vectors::splitAsList(
    as.character(def_alt@unlistData),
    IRanges::togroup(PartitioningByWidth(def_alt))
  ))
  REF <- as.character(def_ref)
  converted <- .geno2geno(NULL, ALT, REF, GT)
  VariantAnnotation::geno(p)$GT <- converted
  callableAlleles(p) <- allele

  return(p)
})

#' Extracts star allele e.g. CYP2D6_12 -> *12
.extractStarAllele <- function(allele) {
  star <- regmatches(allele, regexpr("\\*[0-9]+$", allele))

  return(star)
}

#' Call the phased diplotype for a single sample (column) of the genotype matrix
.callSamplePhasedDiplotype <- function(x, allele, REF, ALT) {
  star <- .extractStarAllele(allele)

  gts <- strsplit(x, "|", fixed = TRUE)
  h1 <- unlist(lapply(gts, function(x) x[1]))
  h2 <- unlist(lapply(gts, function(x) x[2]))
  hm <- rbind(h1, h2)

  call1 <- "*1"
  call2 <- "*1"
  if (all(hm[1, ] == ALT))
      call1 <- star
  if (all(hm[2, ] == ALT)) 
      call2 <- star

  call_string <- paste0(call1, "|", call2)
  
  return(call_string)
}

#' Call the phased diplotypes for all samples in the PGx object
#'
#' Returns a DataFrame of allele calls for every sample in the single allele PGx
#' object.
#' @export
#' @rdname callPhasedDiplotype
#' @return DataFrame with allele calls for every sample in PGx
setMethod("callPhasedDiplotype", "PGx", function(x) {
  if (length(callableAlleles(x)) != 1 & callableAlleles(x) != "") {
    stop("Only a single allele can be called with this method. Run pgxGenotypeCodesToNucleotides(x, allele)")
  }
  allele <- callableAlleles(x)
  allele_gr <- availableHaplotypeRanges(allele, build = pgxBuild(x))
  allele_gr <- allele_gr[GenomicRanges::match(SummarizedExperiment::rowRanges(x), allele_gr), ]
  def_ref <- allele_gr$`Reference Allele`
  def_alt <- allele_gr$`Variant Allele`
  
  # Special case for nested CYP2C19*35 in CYP2C19*2
  if (allele == "CYP2C19*35")
      def_alt <- allele_gr$`Variant Allele2`
  
  # Special case for CYP2C9*8 and CYP2C9*27 -- use ALT defined by sample VCF
  if (allele == "CYP2C9*8" || allele == "CYP2C9*27")
      def_alt <- alt(x)

  gt <- geno(x)$GT
  lst <- asplit(gt, MARGIN = 2)
  calls <- lapply(lst, .callSamplePhasedDiplotype, allele = allele, REF = def_ref, ALT = def_alt)
  DF <- S4Vectors::DataFrame(Call = unlist(calls), row.names = names(calls))
  colnames(DF) <- allele

  return(DF)
})

#' Call phased diplotypes for all samples and callable alleles of a PGx object
#'
#' @export
#' @rdname callPhasedDiplotypes
#' @return DataFrame containing allele calls for all samples of the PGx object
#' for each callable allele
setMethod("callPhasedDiplotypes", "PGx", function(x) {
  x <- getCallableAlleles(x)
  alleles <- callableAlleles(x)
  stopifnot("No callable alleles!" = length(alleles) >= 1 & alleles[1] != "")
  pgx_objs <- lapply(alleles, function(allele) pgxGenotypeCodesToNucleotides(x, allele))
  dfs <- lapply(pgx_objs, callPhasedDiplotype)
  DF <- do.call(cbind, dfs)

  return(DF)
})
