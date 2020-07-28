#' Estimate simATAC simulation parameters
#'
#' Estimate parameters for the simATAC simulation from a real bin by cell
#' input matrix.
#'
#' @param count Either a sparse bin by cell count matrix, or a SingleCellExperiment
#'        object containing count matrix to estimate parameters from.
#' @param object simATACCount object to store estimated parameters and
#'        count matrix in it.
#' @param verbose logical variable. Prints the simulation progress if TRUE.
#'
#' @return simATACCount object containing the estimated counts and parameters.
#'
#' @example
#' count <- getCountFromh5("GSE99172.snap")
#' object <- simATACEstimate(count = count)
#'
#' @seealso
#' \code{\link{simATACEstLibSize}}, \code{\link{simATACEstNonZeroPro}}
#' \code{\link{simATACEstBinMean}}
#' @export
#'
simATACEstimate <- function(count, object = newsimATACCount(), verbose = TRUE) {
  UseMethod("simATACEstimate")
}


#' @rdname simATACEstimate
#' @importFrom methods as
#' @export
#'
simATACEstimate.SingleCellExperiment <- function(count, object = newsimATACCount(), verbose = TRUE) {
  count <- getCountFromSCE(count)
  simATACEstimate.dgCMatrix(as(count, "dgCMatrix"), object, verbose)
}


#' @rdname simATACEstimate
#' @importFrom stats median
#' @importFrom utils write.table
#' @importFrom Matrix colSums
#' @export
#'
simATACEstimate.dgCMatrix <- function(count, object = newsimATACCount(), verbose = TRUE) {

  checkmate::assertClass(object, "simATACCount")

  object <- setParameters(object, nBins=nrow(count))
  count <- count[,which(colSums(count) != 0)]

  lib.size <- colSums(count)
  lib.med <- median(lib.size)
  norm.count <- as(t(t(as.matrix(count)) / lib.size * lib.med), "dgCMatrix")

  if (dim(count)[1] > 0){
    object <- setParameters(object, default = FALSE)
  }

  if (verbose) {message("simATAC is:")}

  if (verbose) {message("...Estimating library size...")}
  object <- simATACEstLibSize(count, object)

  if (verbose) {message("...Estimating non-zero cell proportion...")}
  object <- simATACEstNonZeroPro(norm.count, object)

  if (verbose) {message("...Estimating bin mean...")}
  object <- simATACEstBinMean(norm.count, object)

  return(object)
}


#' Estimate library size parameters
#'
#' A Gaussian mixture distribution is fitted to the log2 transformation
#' of library size values. The estimated parameters are added to the object
#' variable. See \code{\link[mixtools]{normalmixEM}} for details on the fitting.
#'
#' @param count A sparse count matrix to estimate parameters from.
#' @param object simATACCount object to store estimated parameters.
#'
#' @return simATACCount object with updated library size parameter.
#'
#' @importFrom mixtools normalmixEM
#' @importFrom Matrix colSums
#'
simATACEstLibSize <- function(count, object) {

  lib.size <- colSums(count)
  lib.size <- lib.size[lib.size > 0]
  lib.size <- lib.size[!is.na(lib.size)]
  lib.size <- log2(lib.size+1)

  invisible(capture.output(fit <- normalmixEM(lib.size, k = 2, verb = FALSE)))

  lib.mean1 <- unname(fit$mu[1])
  lib.sd1 <- unname(fit$sigma[1])
  lib.mean2 <- unname(fit$mu[2])
  lib.sd2 <- unname(fit$sigma[2])
  lib.prob <- unname(fit$lambda[1])

  object <- setParameters(object,
                          lib.mean1 = lib.mean1,
                          lib.sd1 = lib.sd1,
                          lib.mean2 = lib.mean2,
                          lib.sd2 = lib.sd2,
                          lib.prob = lib.prob)

  return(object)
}


#' Estimate bins' non-zero cell proportion
#'
#' Extract the non-zero cell proportion of each bin among all single cells from
#' the input count matrix.
#'
#' @param count A sparse count matrix to estimate parameters from.
#' @param object simATACCount object to store estimated parameters.
#'
#' @return simATACCount object with updated non-zero cell proportion parameter.
#'
#' @details
#' Vector of non-zero cell proportions of bins is calculated by
#' dividing the number of non-zero entries over the number of all cells
#' for each bin.
#' @importFrom Matrix rowSums
#'
simATACEstNonZeroPro <- function(count, object) {

  non.zero.pro <- rowSums(count != 0)/ncol(count)
  object <- setParameters(object, non.zero.pro = non.zero.pro)

  return(object)
}


#' Estimate bin coverage mean
#'
#' Estimate the coefficients of second degree polynomial relation between bins' non-zero
#' cell proportions and means.
#'
#' @param count A sparse count matrix to estimate parameters from.
#' @param object simATACCount object to store estimated parameters.
#'
#' @return simATACCount object with updated polynomial function parameters.
#'
#' @details
#' Parameters for the second degree polynomial function between bin means and bin
#' non-zero cell proportion are estimated by using \code{\link[stats]{lm}} and
#' \code{\link[stats]{poly}}. simATAC excludes bins having non-zero cell proportion
#' more than 0.8 in modelling the polynomial function.
#'
#' @importFrom stats lm poly coef
#' @importFrom Matrix rowSums
#'
simATACEstBinMean <- function(count, object) {

  non.zero.pro <- simATACget(object, "non.zero.pro")
  bin.mean <- rowSums(count)/ncol(count)

  fit <- lm(bin.mean[which(non.zero.pro < 0.8)]~poly(non.zero.pro[which(non.zero.pro < 0.8)], 2, raw=TRUE))

  c0 <- unname(coef(fit)[1])
  c1 <- unname(coef(fit)[2])
  c2 <- unname(coef(fit)[3])

  object <- setParameters(object, mean.coef0=c0, mean.coef1=c1, mean.coef2=c2)

  return(object)
}
