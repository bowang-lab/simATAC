#' Estimate and simulate simATAC simulation parameters
#'
#' Estimate parameters for the simATAC simulation from a real bin by cell
#' input matrix, and use them for simualting final counts.
#'
#' @param count Either a sparse bin by cell count matrix, or a SingleCellExperiment
#'        object containing count matrix to estimate parameters from.
#' @param object simATACCount object to store estimated parameters and
#'        count matrix in it.
#' @param default Logical variable. Sets default parameters if TRUE.
#' @param verbose logical variable. Prints the simulation progress if TRUE.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{simATACCount} object.
#'
#' @return SingleCellExperiment object containing the estimated counts and parameters.
#'
#'
#' @example
#' count <- getCountFromh5("GSE99172.snap")
#'
#' # simple simulation
#' sim <- simATACGenerate(coutn=t(count))
#'
#'\dontrun{
#' # set nCells parameter
#' sim <- simATACGenerate(count=count, nCells=500)
#'
#' # simulation with default parameters
#' sim <- simATACGenerate(default = TRUE, nCells=500)
#' }
#'
#' @seealso
#' \code{\link{simATACEstimate}}, \code{\link{setParameters}}
#' \code{\link{simATACSimulate}}
#' @export
#'
simATACGenerate <- function(count = NULL, object = newsimATACCount(), default = TRUE, verbose = TRUE, ...) {


  if(is.null(count)){
    default <- FALSE
  }
  if (default == FALSE){
    object <- simATACEstimate(count, object, verbose)
  }

  object <- setParameters(object, ...)
  sim <- simATACSimulate(object)

  return(sim)
}


#' Return count matrix from a SingleCellExperiment object. If count matrix is missing
#' a warning is printed and the first assay is returned.
#'
#' @param sce SingleCellExperiment input object containing counts.
#'
#' @return Sparse matrix containing counts.
#'
#' @example
#' \dontrun{
#' # gets cell by bin sparse matrix
#' count <- getCountFromh5("GSE99172.snap")
#'
#' # create SingleCellExperiment object with bin by cell matrix
#' sce <- SingleCellExperiment(assays = list(counts = t(count)))
#'
#' object <- simATACEstimate(sce)
#' sim <- simATACSimulate(object)
#' }
#' @export
#'
getCountFromSCE <- function(sce) {

  checkmate::assertClass(sce, "SingleCellExperiment")

  if ("counts" %in% SummarizedExperiment::assayNames(sce)) {
    count <- SingleCellExperiment::counts(sce)
  } else {
    warning("counts assay is missing, using the first assay instead")
    count <- SummarizedExperiment::assay(sce)
  }

  return(count)
}


#' Simulate a bin's mean from its non-zero cell proportion
#' from estimated second degree polynomial coefficients.
#'
#' @param x non-zero cell proportion of a bin.
#' @param c0 coefficient of x power 0.
#' @param c1 coefficient of x power 1.
#' @param c2 coefficient of x power 2.
#'
#' @return dependent variable y, which is the simulated bin mean.
#'
simBinMeans <- function(x, c0, c1, c2){

  y=c0+c1*x+c2*x^2
  y = max(0, y)

  return(y)
}


#' Extract count matrix from h5 file
#'
#' @param file Input .h5 (or snap) file to read cell by bin matrix from
#'
#' @return The cell by bin count matrix from .h5 file
#'
#' @example
#' \dontrun{
#' count <- getCountFromh5("GSE99172.snap")
#' }
#'
#' @importFrom Matrix sparseMatrix
#' @importFrom rhdf5 h5read
#' @export
#'
getCountFromh5 <- function(file){

  count.x <- h5read(file, "/AM/5000/idx")
  count.y <- h5read(file, "/AM/5000/idy")
  count.value <- h5read(file, "/AM/5000/count")
  count.barcodeLen <- h5read(file, "/FM/barcodeLen")
  count.binChrom <- h5read(file, "/AM/5000/binChrom")
  count.binStart <- h5read(file, "/AM/5000/binStart")

  nCells <- length(count.barcodeLen)
  nBins <- length(count.binChrom)

  bin.names <- sapply(seq(1, nBins, 1), function(x) paste(count.binChrom[x], count.binStart[x], sep = "_"))
  count <- Matrix::sparseMatrix(count.x, count.y, x = count.value, dims = c(nCells, nBins))

  return(count)
}


#' Convert raw bin by cell matrix in an SingleCellExperiment object into the binary version.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#'
#' @return Sparse matrix containing binary version of the input simulated bin by cell matrix.
#'
#' @example
#' object <- newsimATACCount()
#' sim <- simATACSimulate(object, default = TRUE)
#' count.bin <- simATACgetBinary(sim)
#'
#' @export
#'
simATACgetBinary <- function(sim){

  count.bin <- BiocGenerics::counts(sim)
  count.bin[count.bin > 0] <- 1
  count.bin <- as(count.bin, "dgCMatrix")

  return(count.bin)

}


#' Convert raw bin by cell matrix in an SingleCellExperiment object into the
#' peak by cell matrix by extracting bins having the highest means.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param peak.num Number of peak bins to extract from original bin by cell matrix.
#'
#' @return Sparse matrix containing peak.num bins with the highest bin means.
#'
#' @example
#' object <- newsimATACCount()
#' sim <- simATACSimulate(object, default = TRUE)
#' count.bin <- simATACgetCellByPeak(sim)
#'
#' @importFrom Matrix rowSums
#' @export
#'
simATACgetCellByPeak <- function(sim, peak.num = 5000){

  bin.mean <- rowSums(BiocGenerics::counts(sim))/ ncol(BiocGenerics::counts(sim))
  peak.index <- order(bin.mean, decreasing=TRUE)[1:peak.num]
  count.peak <- as(BiocGenerics::counts(sim)[peak.index,], "dgCMatrix")

  return(count.peak)
}


#' This function returns the chr:start-end name format of the input region from bed file. start
#' and end are coordinated of the begining and end of a specific bin with 5000 length base pairs.
#'
#' @param region input list in the format of [chr, start, end]
#' @param bin.name List of bin names of the raw bin by cell matrix
#'
#' @return List of bin names that has intersection with the input region
#'
getBin <- function(region, bin.name){

  int.div.start <- as.numeric(region[2])%/%5000
  int.div.end <- as.numeric(region[3])%/%5000
  remainder.end <- as.numeric(region[3])%%5000

  start <- int.div.start
  end <- if (remainder.end == 0) int.div.end-1 else int.div.end

  bed.name <- sapply(start:end, function(x)
    paste(region[1], ":", as.character(x*5000+1), "-", as.character((x+1)*5000), sep = ""))

  bed.index <- which(bin.name == bed.name)
  bed.index <- paste(as.character(bed.index), collapse = ':')

  return(bed.index)
}


#' Convert raw bin by cell matrix in an SingleCellExperiment object into the region by cell
#' matrix, with regions defined in the input bed file.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param file.bed Bed file containing chromosome and start and end positions of regions.
#'
#' @return Sparse matrix containing bins having intersection with the regions in the bed file.
#'
#' @example
#' \dontrun{
#' object <- newsimATACCount()
#' sim <- simATACSimulate(object, default = TRUE)
#' count.bin <- simATACgetCellByRegion(sim, file.bed = "file.bed")
#' }
#'
#' @export
#'
simATACgetCellByRegion <- function(sim, file.bed){

  bin.name <- rownames(BiocGenerics::counts(sim))

  bed <- as.data.frame(read.table(file.bed, header = FALSE, sep="", stringsAsFactors=FALSE, quote=""))
  bed <- bed[,1:3]
  colnames(bed) <- c("chr", "start", "end")

  bed.index <- apply(bed, 1, getBin, bin.name = bin.name)
  bed.index <- paste(as.character(bed.index), collapse = ':')

  bed.index <- strsplit(bed.index,":")
  bed.index <- unique(unlist(bed.index))

  count <- BiocGenerics::counts(sim)[as.numeric(bed.index), ]

  return(count)
}
