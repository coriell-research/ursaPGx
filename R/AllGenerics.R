#' @export
setGeneric("pgxGene", function(x) standardGeneric("pgxGene"))

#' @export
setGeneric("pgxBuild", function(x) standardGeneric("pgxBuild"))

#' @export
setGeneric("pgxCallableAlleles", function(x) standardGeneric("pgxCallableAlleles"))

#' @export
setGeneric("pgxCallableAlleles<-", function(x, ..., value) standardGeneric("pgxCallableAlleles<-"))

#' @export
setGeneric("determineCallableAlleles", function(x) standardGeneric("determineCallableAlleles"))

#' @export
setGeneric("pgxReferenceDataFrame", function(x) standardGeneric("pgxReferenceDataFrame"))

#' @export
setGeneric("pgxReferenceDataFrame<-", function(x, ..., value) standardGeneric("pgxReferenceDataFrame<-"))

#' @export
setGeneric("buildReferenceDataFrame", function(x) standardGeneric("buildReferenceDataFrame"))

#' @export
setGeneric("convertGTtoNucleotides", function(x) standardGeneric("convertGTtoNucleotides"))

#' @export
setGeneric("convertGTtoNucleotides<-", function(x, ..., value) standardGeneric("convertGTtoNucleotides<-"))

#' @export
setGeneric("callPhasedDiplotypes", function(x) standardGeneric("callPhasedDiplotypes"))

#' @export
setGeneric("pgxGenotypeMatrix", function(x) standardGeneric("pgxGenotypeMatrix"))

#' @export
setGeneric("detailPhasedCalls", function(x) standardGeneric("detailPhasedCalls"))
