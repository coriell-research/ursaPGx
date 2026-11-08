#' @rdname pgxGene
#' @export
setGeneric("pgxGene", function(x) standardGeneric("pgxGene"))

#' @rdname pgxBuild
#' @export
setGeneric("pgxBuild", function(x) standardGeneric("pgxBuild"))

#' @rdname pgxCallableAlleles
#' @export
setGeneric("pgxCallableAlleles", function(x) standardGeneric("pgxCallableAlleles"))

#' @rdname pgxCallableAlleles
#' @export
setGeneric("pgxCallableAlleles<-", function(x, ..., value) standardGeneric("pgxCallableAlleles<-"))

#' @rdname determineCallableAlleles
#' @export
setGeneric("determineCallableAlleles", function(x) standardGeneric("determineCallableAlleles"))

#' @rdname pgxReferenceDataFrame
#' @export
setGeneric("pgxReferenceDataFrame", function(x) standardGeneric("pgxReferenceDataFrame"))

#' @rdname pgxReferenceDataFrame
#' @export
setGeneric("pgxReferenceDataFrame<-", function(x, ..., value) standardGeneric("pgxReferenceDataFrame<-"))

#' @rdname buildReferenceDataFrame
#' @export
setGeneric("buildReferenceDataFrame", function(x) standardGeneric("buildReferenceDataFrame"))

#' @rdname convertGTtoNucleotides
#' @export
setGeneric("convertGTtoNucleotides", function(x) standardGeneric("convertGTtoNucleotides"))

#' @rdname convertGTtoNucleotides
#' @export
setGeneric("convertGTtoNucleotides<-", function(x, ..., value) standardGeneric("convertGTtoNucleotides<-"))

#' @rdname callPhasedDiplotypes
#' @export
setGeneric("callPhasedDiplotypes", function(x) standardGeneric("callPhasedDiplotypes"))

#' @rdname pgxGenotypeMatrix
#' @export
setGeneric("pgxGenotypeMatrix", function(x) standardGeneric("pgxGenotypeMatrix"))

#' @rdname detailPhasedCalls
#' @export
setGeneric("detailPhasedCalls", function(x) standardGeneric("detailPhasedCalls"))
