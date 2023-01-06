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

#' Extract a list of VRanges for all defined haplotypes of the PGx object
#' 
#' Return a list of VRanges objects for all of the defined haplotype ranges of
#' the given PGx object.
#' @export
#' @rdname extractHaplotypeRanges
#' @return list of VRanges object containing all allele definitions for the 
#' given PGx object
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
    vrl <- extractHaplotypeRanges(x)
    q <- SummarizedExperiment::rowRanges(x)
    callable <- lapply(vrl, function(x) .isCallable(q, s = x))
    callable <- Filter(Negate(is.null), callable)
    stopifnot("There are no callable alleles" = length(callable) >= 1)
    callableAlleles(x) <- names(callable)
    
    return(x)
})
