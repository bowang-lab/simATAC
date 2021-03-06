#' simATAC simulation
#'
#' Simulate bin by cell count matrix from a sparse single-cell ATAC-seq bin by cell
#' input using simATAC methods.
#'
#' @param object simATACCount object with simulation parameters. See \code{\link{simATACCount}} 
#'        for details.
#' @param verbose Logical variable. Prints the simulation progress if TRUE.
#' @param ... Any additional parameter settings to override what is provided in
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
#' simATAC provides the option to manually adjust each of the \code{simATACCount}
#' object parameters by calling \code{\link{setParameters}}. See examples for a demonstration
#' of how this can be used.
#'
#' The simulation involves the following steps:
#' \enumerate{
#'     \item Set up simulation parameters
#'     \item Set default parameters if needed
#'     \item Set up SingleCellExperiment object
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
#'             \item{BinMean}{The simulated bins' means.}
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
#' Code: \url{https://github.com/bowang-lab/simATAC}
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
  if (verbose) {message("...updating parameters...")}
  object <- setParameters(object, ...)
  
  validObject(object)
  
  default <- simATACget(object, "default")
  if (default == TRUE){
    if (verbose) {message("...setting default parameters...")}
    object <- setDefautParameters(object)
  }
  
  # Set random seed in each run
  seed <- simATACget(object, "seed")
  set.seed(seed)
  
  nCells <- simATACget(object, "nCells")
  nBins <- simATACget(object, "nBins")
  
  if (verbose) {message("...setting up SingleCellExperiment object...")}
  
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
  
  if (verbose) {message("...simulating library size...")}
  sim <- simATACSimLibSize(object, sim)
  gc()
  
  if (verbose) {message("...simulating non-zero cell proportion...")}
  object <- simATACSimZeroEntry(object)
  gc()
  
  if (verbose) {message("...simulating bin mean...")}
  sim <- simATACSimBinMean(object, sim)
  gc()
  
  if (verbose) {message("...generating final counts...")}
  sim <- simATACSimTrueCount(object, sim)
  
  colnames(BiocGenerics::counts(sim)) <- cell.names
  rownames(BiocGenerics::counts(sim)) <- bin.names
  
  if (verbose) {message("...Done...")}
  
  return(sim)
}


#' Calculate the polynomial coefficients from an existing dataset GSE99172,
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
  
  path <- system.file("extdata", "GSE99172.txt", package="simATAC")
  data <- read.table(path, sep=",")
  colnames(data) <- c("bin.mean", "non.zero.pro")
  
  fit <- lm(as.double(as.numeric(data[which(data$non.zero.pro < 0.8),1]))~poly(as.double(as.numeric(data[which(data$non.zero.pro < 0.8),2])), 2, raw=TRUE))
  c0 <- unname(coef(fit)[1])
  c1 <- unname(coef(fit)[2])
  c2 <- unname(coef(fit)[3])
  
  object <- setParameters(object,
                          non.zero.pro=as.numeric(data[,2]),
                          mean.coef0=c0,
                          mean.coef1=c1,
                          mean.coef2=c2)
  
  return(object)
}


#' Set bin names based on the species variable of the simATACCount object.
#'
#' @param object simATACCount object with simulation parameters.
#'
#' @return List of bin names.
#'
getBinNames <- function(object){
  
  species <- simATACget(object, "species")
  nBins <- simATACget(object, "nBins")
  bin.coordinate.file <- simATACget(object, "bin.coordinate.file")
  ref <- c('hg19', 'hg38', 'mm9', 'mm10')
  file <- ""
  
  if  (!(species %in% ref)){
    message(species, " referene is not supported by current version of simATAC. simATAC supports hg19, hg38, mm9, and mm10 references. Please give a file of bin information consistent with your input data with three columns and header of \"chr start end\" as the bin.coordinate.file parameter. If you don't have a file containing the information of bins, simATAC considers the bin.coordinate.file parameter as \"None\". In this case, you will not be able to get the coordinate information of bins. Please make sure the \"species\" parameter of the simATACCount object is set correctly.")
    bin.names <- paste0("Bin", seq_len(nBins))
  }
  
  if (bin.coordinate.file != 'None'){
    file <- bin.coordinate.file
    data <- read.table(file, header = TRUE)
    bin.names <- paste(tolower(data$chr), ":", data$start, "-", data$end, sep = "")
  }else if (species %in% ref){
    file <- system.file("extdata", paste(species, "_genome_coordinates.txt", sep=""), package = "simATAC")
    data <- read.table(file, header = TRUE)
    bin.names <- paste(tolower(data$chr), ":", data$start, "-", data$end, sep = "")
  }
  
  if (length(bin.names) != nBins){
    message("Your data has different number of bins compared to the provided genome positions. Please give a file of bin information consistent with your input data with three columns and header of \"chr start end\" as the bin.coordinate.file parameter. If you don't give a file containing the information of bins, simATAC considers the bin.coordinate.file parameter as \"None\" and names the bins {Bin1 to BinX} with X number of bins. In this case, you wont be able to get the coordinate information of bins. Please make sure the \"species\" parameter of the simATACCount object is set correctly.")
    bin.names <- paste0("Bin", seq_len(nBins))
  }
  
  
  
  
  
  return(bin.names)
}


#' Simulate library sizes by sampling from a bimodal Gaussian mixture model based
#' on the estimated mus and sigmas and weight of its components.
#'
#' @param object simATACCount object with simulation parameters.
#' @param sim SingleCellExperiment object containing simulation parameters.
#'
#' @return SingleCellExperiment object with updated simulation library sizes.
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
#' The extracted non-zero cell proportion from the original bin by cell 
#' matrix is used as probability to simulate cells with zero counts, using 
#' Bernoulli distribution.
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


#' Simulate bin means from estimated polynomial coefficients using
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


#' Adjust the simulated library sizes and bin means to generate final counts,
#' using Poisson distribution, with optional Gaussian noise.
#'
#' @param object simATACCount object with simulation parameters.
#' @param sim SingleCellExperiment object containing simulation parameters.
#'
#' @return SingleCellExperiment with simulated final counts.
#'
#' @importFrom SummarizedExperiment rowData colData<- colData colData<- assays assays<-
#' @importFrom stats rpois rnorm
#'
simATACSimTrueCount <- function(object, sim) {
  
  nCells <- simATACget(object, "nCells")
  nBins <- simATACget(object, "nBins")
  noise.mean <- simATACget(object, "noise.mean")
  noise.sd <- simATACget(object, "noise.sd")
  sparse.fac <- simATACget(object, "sparse.fac")
  bin.mean <- rowData(sim)$BinMean
  lib.size <- colData(sim)$LibSize
  mean.sum <- sum(bin.mean)
  multi.fac <- bin.mean/mean.sum
  
  counts <- as(sapply(1:nCells, 
                      function(i)
                        rpois(n=nBins, lib.size[i]*multi.fac*sparse.fac)),
               "dgCMatrix")
  gc()
  
  if (noise.sd > 0){
    
    end <- round(nCells/3)
    counts1 <- as(sapply(1:end,
                         function(i)
                           round(counts[,i]+rnorm(n=nBins, mean = noise.mean, sd = noise.sd))),
                  "dgCMatrix")
    counts1[counts1 <= 0] <- 0
    gc()
    
    start <- end+1
    end <- round(2*nCells/3)
    counts2 <- as(sapply(start:end,
                         function(i)
                           round(counts[,i]+rnorm(n=nBins, mean = noise.mean, sd = noise.sd))),
                  "dgCMatrix")
    counts2[counts2 <= 0] <- 0
    gc()
    
    start <- end+1
    counts3 <- as(sapply(start:nCells,
                         function(i)
                           round(counts[,i]+rnorm(n=nBins, mean = noise.mean, sd = noise.sd))),
                  "dgCMatrix")
    counts3[counts3 <= 0] <- 0
    gc()
    
    counts <- cbind(counts1, counts2, counts3)
  }
  
  assays(sim, withDimnames = FALSE)$counts <- counts
  
  return(sim)
}
