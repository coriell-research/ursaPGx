% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PGx-methods.R
\name{determineCallableAlleles,PGx-method}
\alias{determineCallableAlleles,PGx-method}
\title{Return a vector of allele names that are able to be called for the given PGx object}
\usage{
\S4method{determineCallableAlleles}{PGx}(x)
}
\value{
PGx object with pgxCallableAlleles slot filled
}
\description{
Callable alleles are alleles where all defined positions are present in the
sample VCF. The function will also add the vector of the callable alleles to
the pgxCallableAlleles slot of the PGx object
}
