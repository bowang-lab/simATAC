% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ExtraFunctions.R
\name{simATACgetFeatureByCell}
\alias{simATACgetFeatureByCell}
\title{Convert raw bin by cell matrix in an SingleCellExperiment object into the region by cell
matrix, with regions defined in the input BED file.}
\usage{
simATACgetFeatureByCell(sim, file.bed)
}
\arguments{
\item{sim}{SingleCellExperiment object containing simulation parameters.}

\item{file.bed}{BED file containing chromosome, start, and end positions of regions.}
}
\value{
Sparse matrix containing bins having intersection with the regions in the BED file.
}
\description{
Convert raw bin by cell matrix in an SingleCellExperiment object into the region by cell
matrix, with regions defined in the input BED file.
}
