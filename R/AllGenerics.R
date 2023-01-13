#' @export
setGeneric("pgxGene", function(x) standardGeneric("pgxGene"))

#' @export
setGeneric("pgxBuild", function(x) standardGeneric("pgxBuild"))

#' @export
setGeneric("callableAlleles", function(x) standardGeneric("callableAlleles"))

#' @export
setGeneric("callableAlleles<-", function(x, ..., value) standardGeneric("callableAlleles<-"))

#' @export
setGeneric("getCallableAlleles", function(x) standardGeneric("getCallableAlleles"))

#' @export
setGeneric("extractHaplotypeRanges", function(x) standardGeneric("extractHaplotypeRanges"))

#' @export
setGeneric("pgxGenotypeCodesToNucleotides", function(x, allele, ...) standardGeneric("pgxGenotypeCodesToNucleotides"))

#' @export
setGeneric("callPhasedDiplotype", function(x) standardGeneric("callPhasedDiplotype"))

#' @export
setGeneric("callPhasedDiplotypes", function(x) standardGeneric("callPhasedDiplotypes"))