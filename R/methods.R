#' Show method for PGx objects
#' @rdname PGx-class
#' @export
setMethod("show", "PGx", function(object) {
  cat(
    "Class: PGx\n",
    "PGx gene: ", object@pgxGene, " \n",
    "PGx build: ", object@pgxBuild, " \n",
    "Number of samples: ", ncol(object), " \n",
    sep = ""
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
