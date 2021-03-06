% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtraFunctions.R
\name{simATACGenerate}
\alias{simATACGenerate}
\title{Estimate and simulate simATAC simulation parameters}
\usage{
simATACGenerate(
  count = NULL,
  object = newsimATACCount(),
  default = TRUE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{count}{Either a sparse bin by cell count matrix, or a SingleCellExperiment
object containing the count matrix to estimate parameters from.}

\item{object}{A simATACCount object to store estimated parameters and
the count matrix in it.}

\item{default}{Logical variable. Sets default parameters if TRUE.}

\item{verbose}{Logical variable. Prints the simulation progress if TRUE.}

\item{...}{Any additional parameter settings to override what is provided in
\code{simATACCount} object.}
}
\value{
SingleCellExperiment object containing the estimated counts and parameters.
}
\description{
Estimate parameters for the simATAC simulation from a real bin by cell
input matrix, and use them for simualting final counts.
}
\seealso{
\code{\link{simATACEstimate}}, \code{\link{setParameters}}
\code{\link{simATACSimulate}}
}
