#' PGx class object
#'
#' The PGx class inherits from \code{\link[VariantAnnotation]{VCF}}. Like the
#' \code{CollapsedVCF} class it is expected that users will not create instances
#' of the PGx class but instead the PGx class will be created by one of the
#' constructor functions, i.e. \code{readPGx}.
#'
#' @slot pgxGene The PGx gene used to construct the PGx object. See \code{pgxGenes()} for available genes.
#' @slot pgxBuild The genome build of the PGx object.
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
    pgxReferenceDataframe = "data.frame"
  ),
  contains = "CollapsedVCF"
)

#' @export
#' @importFrom GenomicRanges GRanges
#' @importClassesFrom VariantAnnotation VCF
PGx <- function(vcf = VCF(collapsed = TRUE), pgxBuild = "", pgxGene = "",
                pgxCallableAlleles = "", pgxReferenceDataframe = data.frame()) {
  .PGx(vcf, pgxGene = pgxGene, pgxBuild = pgxBuild, pgxCallableAlleles = pgxCallableAlleles,
       pgxReferenceDataframe = pgxReferenceDataframe)
}

setValidity("PGx", function(object) {
    if (!is(object@pgxGene)[[1]] == "character") {
        return("pgxGene must be non-empty")
        }
    if (!is(object@pgxBuild)[[1]] == "character") {
        return("pgxBuild must be non-empty")
    }
    if (length(object@pgxCallableAlleles) == 1) {
        if (object@pgxCallableAlleles != "" && !object@pgxCallableAlleles %in% availableHaplotypes()) {
            return("callableAlleles must be empty or one of availableHaplotypes()")
        }
    }
    if (length(object@pgxCallableAlleles) > 1) {
        if (!all(object@pgxCallableAlleles %in% availableHaplotypes())) {
            return("pgxCallableAlleles must be in availableHaplotypes()")
        }
    }
    TRUE
})
