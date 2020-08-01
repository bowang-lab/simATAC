#' simATAC simulation
#'
#' Simulate bin by cell count matrix from a sparse single-cell ATAC-seq bin by cell
#' input using simATAC methods.
#'
#' @param object simATACCount object with simulation parameters.
#'        See \code{\link{simATACCount}} for details.
#' @param verbose logical variable. Prints the simulation progress if TRUE.
#' @param ... any additional parameter settings to override what is provided in
#'        \code{simATACCount} object.
#'
#' @return SingleCellExperiment object containing the simulated counts.
#'
#' @example
#' # simple simulation
#' count <- getCountFromh5("GSE99172.snap")
#' object <- simATACEstimate(count = count)
#' sim <- simATACSimulate(object)
#'
#' object <- simATACEstimate(count = count)
#' object <- setParameters(object, nCells=500)
#' sim <- simATACSimulate(object)
#'
#' object <- simATACEstimate(count = count)
#' object <- setParameters(object)
#' sim <- simATACSimulate(object, nCells=500, noise.mean=-0.3, noise.sd=0.3)
#'
#' @details
#'simATAC provides the option to manually adjust each of the \code{simATACCount}
#' object parameter by calling \code{\link{setParameters}}. See examples for a demonstration
#' of how this can be used.
#'
#' The simulation involves the following steps:
#' \enumerate{
#'     \item Set up simulation object
#'     \item Simulate library sizes
#'     \item Simulate non zero entries
#'     \item Simulate bin means
#'     \item Create final synthetic counts
#' }
#'
#' The final output is a \code{\link[SingleCellExperiment]{SingleCellExperiment}} object that
#' contains the simulated count matrix. The parameters are stored in the \code{\link{colData}}
#' (for cell specific information), \code{\link{rowData}} (for bin specific information) or
#' \code{\link{assays}} (for bin by cell matrix) slots. This additional information includes:
#' \describe{
#'     \item{\code{rowData}}{
#'         \describe{
#'             \item{BinMean}{The simulated bin means.}
#'         }
#'     }
#'     \item{\code{colData}}{
#'         \describe{
#'             \item{LibSize}{The simulated library size values.}
#'         }
#'     }
#'     \item{\code{assays}}{
#'         \describe{
#'             \item{counts}{The sparse matrix containing simulated counts.}
#'         }
#'     }
#' }
#'
#' @references
#' Navidi Ghaziani Z. (2020).
#'
#' Paper: \url{}
#'
#' Code: \url{}
#'
#' @seealso
#' \code{\link{simATACSimLibSize}}, \code{\link{simATACSimZeroEntry}},
#' \code{\link{simATACSimBinMean}}, \code{\link{simATACSimTrueCount}}
#'
#' @importFrom SummarizedExperiment rowData colData colData<- assays
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom utils capture.output
#' @export
#'
simATACSimulate <- function(object = newsimATACCount(),
                            verbose = TRUE,
                            ...) {
  checkmate::assertClass(object, "simATACCount")

  if (verbose) {message("simATAC is:")}
  if (verbose) {message("...Updating parameters...")}
  object <- setParameters(object, ...)

  validObject(object)

  default <- simATACget(object, "default")
  if (default == TRUE){
    if (verbose) {message("...Setting default parameters...")}
    object <- setDefautParameters(object)
  }

  # Set random seed in each run
  seed <- simATACget(object, "seed")
  set.seed(seed)

  nCells <- simATACget(object, "nCells")
  nBins <- simATACget(object, "nBins")

  if (verbose) {message("...Setting up SingleCellExperiment object...")}

  cell.names <- paste0("Cell", seq_len(nCells))
  bin.names <- getBinNames(object)

  # Create SingleCellExperiment object to store simulation data
  cells <- data.frame(Cell = cell.names)
  rownames(cells) <- cell.names
  bins <- data.frame(Bin = bin.names)
  rownames(bins) <- bin.names

  sim <- SingleCellExperiment(rowData = bins,
                              colData = cells,
                              metadata = list(Params = object))

  if (verbose) {message("...Simulating library size...")}
  sim <- simATACSimLibSize(object, sim)
  gc()

  if (verbose) {message("...Simulating non-zero cell proportion...")}
  object <- simATACSimZeroEntry(object)
  gc()

  if (verbose) {message("...Simulating bin mean...")}
  sim <- simATACSimBinMean(object, sim)
  gc()

  if (verbose) {message("...Generating final counts...")}
  sim <- simATACSimTrueCount(object, sim)

  colnames(BiocGenerics::counts(sim)) <- cell.names
  rownames(BiocGenerics::counts(sim)) <- bin.names

  if (verbose) {message("...Done...")}

  return(sim)
}


#' Calculate the binomial coefficients from an existing dataset GSE99172,
#' for default values.
#'
#' @param object simATACCount object with simulation parameters.
#'
#' @return simATACCount object variable with default variables.
#'
#' @importFrom stats lm poly
#' @importFrom utils read.table
#'
setDefautParameters <- function(object){

  path <- system.file("extdata", "GSE99172.txt", package = "simATAC")
  data <- read.table(path, sep=",")
  colnames(data) <- c("bin.mean", "non.zero.pro")

  fit <- lm(as.double(as.numeric(data[which(data$non.zero.pro < 0.8),1]))~poly(as.double(as.numeric(data[which(data$non.zero.pro < 0.8),2])), 2, raw=TRUE))
  c0 <- unname(coef(fit)[1])
  c1 <- unname(coef(fit)[2])
  c2 <- unname(coef(fit)[3])

  object <- setParameters(object,
                        non.zero.pro = as.numeric(data[,2]),
                        mean.coef0=c0,
                        mean.coef1=c1,
                        mean.coef2=c2)

  return(object)
}


#' Set bin names based on the species parameter of the simATACCount object.
#'
#' @param object simATACCount object with simulation parameters.
#'
#' @return List of bin names.
#'
getBinNames <- function(object){

  species <- simATACget(object, "species")
  file = ""

  if(species == "human"){
    file <- system.file("extdata", "Human_genome_coordinates.txt", package = "simATAC")
  }else if(species == "mouse"){
    file <- system.file("extdata", "Mouse_genome_coordinates.txt", package = "simATAC")
  }

  data <- read.table(file)
  bin.names <- paste(tolower(data$chr), ":", data$start, "-", data$end, sep = "")

  return(bin.names)
}


#' Simulate library sizes by sampling from a bimodal Gaussian mixture model based
#' on the estimated mus and sigmas and weight of its components.
#'
#' @param object simATACCount object with simulation parameters.
#' @param sim SingleCellExperiment object containing simulation parameters.
#'
#' @return SingleCellExperiment object with updated simulation library size.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats rnorm
#'
simATACSimLibSize <- function(object, sim) {

  nCells <- simATACget(object, "nCells")
  lib.mean1 <- simATACget(object, "lib.mean1")
  lib.mean2 <- simATACget(object, "lib.mean2")
  lib.sd1 <- simATACget(object, "lib.sd1")
  lib.sd2 <- simATACget(object, "lib.sd2")
  lib.prob <- simATACget(object, "lib.prob")

  components <- sample(1:2, prob=c(lib.prob, 1-lib.prob), size=nCells, replace=TRUE)
  mus <- c(lib.mean1, lib.mean2)
  sds <- c(lib.sd1, lib.sd2)
  lib.size <- rnorm(n=nCells, mean=mus[components], sd=sds[components])

  lib.size <- 2^lib.size-1
  colData(sim)$LibSize <- lib.size

  return(sim)
}


#' Simulate zero and non-zero entries of simATAC simulation count matrix.
#' The extracted non-zero proportion from original bin by cell matrix is
#' used as probability to set cells to zero using Bernoulli distribution.
#'
#' @param object simATACCount object with simulation parameters.
#'
#' @return simATACCount object with updated simulation parameters.
#'
#' @importFrom SummarizedExperiment rowData colData
#' @importFrom stats rbinom
#' @importFrom Matrix Matrix
#'
simATACSimZeroEntry <- function(object) {

  nCells <- simATACget(object, "nCells")
  nBins <- simATACget(object, "nBins")
  non.zero.pro <- simATACget(object, "non.zero.pro")

  object <- setParameters(object, non.zero.pro = rbinom(nBins, nCells, as.numeric(non.zero.pro))/nCells)

  return(object)
}


#' Simulates bin means from estimated polynomial coefficients using
#' \code{\link{simBinMeans}}.
#'
#' @param object simATACCount object with simulation parameters.
#' @param sim SingleCellExperiment object containing simulation parameters.
#'
#' @return SingleCellExperiment object with updated simulated bin means.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
simATACSimBinMean <- function(object, sim) {

  sim.non.zero.pro <- simATACget(object, "non.zero.pro")

  c0 <- simATACget(object, "mean.coef0")
  c1 <- simATACget(object, "mean.coef1")
  c2 <- simATACget(object, "mean.coef2")

  rowData(sim)$BinMean <- sapply(sim.non.zero.pro, simBinMeans, c0=c0, c1=c1, c2=c2)

  return(sim)
}


#' Adjust the simulated library size and bin means to generate final counts
#' using Poisson distribution, with optional Gaussian noise.
#'
#' @param object simATACCount object with simulation parameters.
#' @param sim SingleCellExperiment object containing simulation parameters.
#'
#' @return SingleCellExperiment with simulated final counts.
#'
#' @importFrom SummarizedExperiment rowData colData<- colData colData<-
#' @importFrom stats rpois rnorm
#'
simATACSimTrueCount <- function(object, sim) {

  nCells <- simATACget(object, "nCells")
  nBins <- simATACget(object, "nBins")
  noise.mean <- simATACget(object, "noise.mean")
  noise.sd <- simATACget(object, "noise.sd")
  bin.mean <- rowData(sim)$BinMean
  lib.size <- colData(sim)$LibSize
  mean.sum <- sum(bin.mean)
  multi.fac <- bin.mean/mean.sum

  BiocGenerics::counts(sim) <- as(sapply(1:nCells, 
                                         function(i)
                                           rpois(n=nBins, lib.size[i]*multi.fac)),
                                  "dgCMatrix")

  if (noise.sd > 0){
    BiocGenerics::counts(sim) <-as(round(BiocGenerics::counts(sim) +
                                           rnorm(nBins*nCells, mean = noise.mean, sd = noise.sd)),
                                   "dgCMatrix")
    BiocGenerics::counts(sim)[BiocGenerics::counts(sim) <= 0] <- 0
  }

  return(sim)
}
