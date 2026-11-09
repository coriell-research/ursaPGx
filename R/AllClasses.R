#' PGx class object
#'
#' The PGx class inherits from \code{\link[VariantAnnotation]{VCF}}. Like the
#' \code{CollapsedVCF} class it is expected that users will not create instances
#' of the PGx class but instead the PGx class will be created by one of the
#' constructor functions, i.e. \code{readPGx()}.
#'
#' @slot pgxGene The PGx gene used to construct the PGx object. See \code{availableGenes()} for available genes.
#' @slot pgxBuild The genome build of the PGx object.
#' @slot pgxCallableAlleles Character vector of all star alleles that are able to be called
#' @slot pgxReferenceDataFrame DataFrame containing reference bases at each position for all callable alleles
#' @slot hasCallableAlleles Logical. Has \code{determineCallableAlleles()} been run?
#' @slot hasReferenceDataFrame Logical. Has \code{buildReferenceDataFrame()} been run?
#' @param vcf A \code{CollapsedVCF} object.
#' @param pgxBuild The genome build. One of "GRCh38" or "GRCh37".
#' @param pgxGene The PGx gene. See \code{availableGenes()}.
#' @param pgxCallableAlleles Character vector of callable star alleles.
#' @param pgxReferenceDataFrame DataFrame of reference haplotype definitions.
#' @param hasCallableAlleles Logical. Have callable alleles been determined?
#' @param hasReferenceDataFrame Logical. Has the reference DataFrame been built?
#' @param object A PGx object.
#' @rdname PGx
#' @aliases PGx-class
#' @export
#' @import methods
#' @importClassesFrom VariantAnnotation CollapsedVCF
.PGx <- setClass(
  "PGx",
  slots = representation(
    pgxGene = "character",
    pgxBuild = "character",
    pgxCallableAlleles = "character",
    pgxReferenceDataFrame = "DFrame",
    hasCallableAlleles = "logical",
    hasReferenceDataFrame = "logical"
  ),
  contains = "CollapsedVCF"
)

#' @rdname PGx
#' @export
#' @importFrom GenomicRanges GRanges
#' @importFrom VariantAnnotation VCF
#' @importClassesFrom VariantAnnotation VCF
PGx <- function(
        vcf = VCF(collapsed = TRUE), 
        pgxBuild = "", 
        pgxGene = "",
        pgxCallableAlleles = "", 
        pgxReferenceDataFrame = S4Vectors::DataFrame(),
        hasCallableAlleles = FALSE,
        hasReferenceDataFrame = FALSE) {
    
  .PGx(vcf, pgxGene = pgxGene, pgxBuild = pgxBuild, 
       pgxCallableAlleles = pgxCallableAlleles,
       pgxReferenceDataFrame = pgxReferenceDataFrame,
       hasCallableAlleles = hasCallableAlleles, 
       hasReferenceDataFrame = hasReferenceDataFrame)
}

setValidity("PGx", function(object) {
    if (!is(object@pgxGene)[[1]] == "character") {
        return("pgxGene must be non-empty")
        }
    if (!is(object@pgxBuild)[[1]] == "character") {
        return("pgxBuild must be non-empty")
    }
    build <- if (length(object@pgxBuild) == 1 && object@pgxBuild %in% c("GRCh38", "GRCh37")) {
        object@pgxBuild
    } else {
        "GRCh38"
    }
    if (length(object@pgxCallableAlleles) == 1) {
        if (object@pgxCallableAlleles != "" && !object@pgxCallableAlleles %in% availableHaplotypes(build)) {
            return("callableAlleles must be empty or one of availableHaplotypes()")
        }
    }
    if (length(object@pgxCallableAlleles) > 1) {
        if (!all(object@pgxCallableAlleles %in% availableHaplotypes(build))) {
            return("pgxCallableAlleles must be in availableHaplotypes()")
        }
    }
    TRUE
})
