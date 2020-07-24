#' The simATACCount class
#'
#' S4 class that holds parameters for the count matrix of simATAC simulation.
#'
#' @section Parameters:
#'
#' simATAC simulation parameters:
#'
#' \describe{
#'     \item{\code{nBins}}{The bin number to simulate.}
#'     \item{\code{nCells}}{The cell number to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{[default]}}{The logical variable whether to use default parameters
#'      (TRUE) or learn from data (FALSE)}
#'      \item{\code{[species]}}{The string indicating the species of the input cells}
#'     \item{\emph{Library size parameters}}{
#'         \describe{
#'             \item{\code{lib.mean1}}{Mean parameter for the first
#'             component of library size bimodal Gaussian distribution.}
#'             \item{\code{lib.mean2}}{Mean parameter for the second
#'             component of library size bimodal Gaussian distribution}
#'             \item{\code{lib.sd1}}{Standard deviation parameter for the
#'             first component of library size bimodal Gaussian distribution.}
#'             \item{\code{lib.sd2}}{Standard deviation parameter for the
#'             second component of library size bimodal Gaussian distribution.}
#'             \item{\code{lib.prob}}{Probability parameter for the first component
#'             in bimodal Gaussian distribution. The probability for the second
#'             component is 1-lib.prob.}
#'         }
#'     }
#'     \item{\emph{Zero entry parameters}}{
#'         \describe{
#'             \item{\code{non.zero.pro}}{The proportion of non-zero cells per
#'             bin in the original count matrix}
#'             \item{\code{mean.coef0}}{Estimated coefficient of power zero
#'             variable in the polynomial function}
#'             \item{\code{mean.coef1}}{Estimated coefficient of power one
#'             variable in the polynomial function}
#'             \item{\code{mean.coef2}}{Estimated coefficient of power two
#'             variable in the polynomial function}
#'         }
#'     }
#'     \item{\emph{[noise]}}{
#'         \describe{
#'             \item{\code{[noise.mean]}}{Gaussian mean to be added as noise to the final simulated
#'    counts}
#'             \item{\code{[noise.sd]}}{Gaussian standard deviation to be added as noise to the final simulated
#'    counts}
#'         }
#'     }
#' }
#'
#' The parameters not shown in brackets can be estimated from real data using
#' \code{\link{simATACEstimate}}. For details of the simATAC simulation
#' see \code{\link{simATACSimulate}}.
#'
#' @name simATACCount
#' @rdname simATACCount
#' @exportClass simATACCount
#' @aliases simATACCount-class
setClass("simATACCount",
         slots = c(nBins = "numeric",
                   nCells = "numeric",
                   seed = "numeric",
                   default = "logical",
                   species = "character",
                   lib.mean1 = "numeric",
                   lib.mean2 = "numeric",
                   lib.sd1 = "numeric",
                   lib.sd2 = "numeric",
                   lib.prob = "numeric",
                   non.zero.pro = "numeric",
                   mean.coef0 = "numeric",
                   mean.coef1 = "numeric",
                   mean.coef2 = "numeric",
                   noise.mean = "numeric",
                   noise.sd = "numeric"),
         prototype = prototype(nBins = 642098,
                               nCells = 500,
                               seed = sample(seq_len(1e6), 1),
                               default = TRUE,
                               species = "human",
                               lib.mean1 = 13.60503,
                               lib.mean2 = 14.93826,
                               lib.sd1 = 1.745264,
                               lib.sd2 = 1.009923,
                               lib.prob = 0.5257138,
                               non.zero.pro = 1,
                               mean.coef0 = 0.002822035,
                               mean.coef1 = 0.6218985,
                               mean.coef2 = 1.976122,
                               noise.mean = 0,
                               noise.sd = 0))
