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
setGeneric("pgxReferenceDataframe", function(x) standardGeneric("pgxReferenceDataframe"))

#' @export
setGeneric("pgxReferenceDataframe<-", function(x, ..., value) standardGeneric("pgxReferenceDataframe<-"))

#' @export
setGeneric("buildReferenceDataframe", function(x) standardGeneric("buildReferenceDataframe"))

#' @export
setGeneric("convertGTtoNucleotides", function(x) standardGeneric("convertGTtoNucleotides"))

#' @export
setGeneric("convertGTtoNucleotides<-", function(x, ..., value) standardGeneric("convertGTtoNucleotides<-"))

#' @export
setGeneric("callPhasedDiplotypes", function(x) standardGeneric("callPhasedDiplotypes"))
