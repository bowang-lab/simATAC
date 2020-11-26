#' Create a new simATACCount object.
#'
#' @param ... Variables to set simATACCount object parameters.
#'
#' @return The new object from class simATACCount.
#'
#' @example
#' object <- newsimATACCount()
#'
#' @importFrom methods new
#' @title newsimATACCount
#' @name newsimATACCount
#' @export
#'
newsimATACCount <- function(...) {

  object <- new("simATACCount")
  object <- setParameters(object, ...)

  return(object)
}


#' Set parameters
#'
#' Set input parameters of the input simATACCount object.
#'
#' @param object Input simATACCount object.
#' @param update Parameter name and value to update in the object variable.
#' @param ... Parameters to set to object variable.
#'
#' @return simATACCount object with updated parameters.
#'
#' @example
#' object <- newsimATACCount()
#' object <- setParameters(object, nCells = 500, noise.mean = -0.3, noise.sd = 0.3)
#'
#'
#' @importFrom methods slot<- validObject
#' @export
#'
setParameters <- function(object, update = NULL, ...){

  checkmate::assertClass(object, classes = "simATACCount")
  checkmate::assertList(update, null.ok = TRUE)

  update <- c(update, list(...))

  if (length(update) > 0) {
    for (name in names(update)) {
      value <- update[[name]]
      checkmate::assertString(name)

      slot(object, name) <- value

      validObject(object)
    }
  }

  return(object)
}


#' Get simATACCount parameters
#'
#' Get the value of input variables from simATACCount object.
#'
#' @param object Input simATACCount object.
#' @param names List of parameter names.
#'
#' @return List of parameter values.
#'
#' @example
#' object <- newsimATACCount()
#' params <- getParameters(object, c(nCells, species, nBins))
#'
#' @export
#'
getParameters <- function(object, names) {

  checkmate::assertClass(object, classes = "simATACCount")
  checkmate::assertCharacter(names, min.len = 1, any.missing = FALSE)

  parameter.list <- lapply(names, simATACget, object = object)
  names(parameter.list) <- names

  return(parameter.list)
}


#' Get a simATACCount parameter
#'
#' Get the value of a single variable from input simATACCount object.
#'
#' @param object Input simATACCount object.
#' @param name Name of the parameter.
#'
#' @return Value of the input parameter.
#'
#' @example
#' object <- newsimATACCount()
#' nBins <- simATACget(object, "nBins")
#'
#' @importFrom methods slot
#' @export
#'
simATACget <- function(object, name) {
  slot(object, name)
}


#' Check the validity of parameters
#'
#' Check the validity of simATACCount object parameters.
#'
#' @param object Input simATACCount object.
#'
#' @return TRUE if all parameters are valid, otherwise FALSE with the list of invalid
#' parameters following.
#'
#' @importFrom checkmate checkInt checkNumber checkNumeric checkFlag checkCharacter checkLogical
#'
setValidity("simATACCount", function(object) {

  name.list <- getParameters(object, c(slotNames(object)))

  checks <- c(nBins = checkInt(name.list$nBins, lower = 1),
              nCells = checkInt(name.list$nCells, lower = 1),
              species = checkCharacter(name.list$species, min.len = 1),
              seed = checkInt(name.list$seed),
              default = checkLogical(name.list$default),
              lib.mean1 = checkNumber(name.list$lib.mean1, lower = 0),
              lib.mean2 = checkNumber(name.list$lib.mean2, lower = 0),
              lib.sd1 = checkNumber(name.list$lib.sd1, lower = 0),
              lib.sd2 = checkNumber(name.list$lib.sd2, lower = 0),
              lib.prob = checkNumber(name.list$lib.prob, lower = 0, upper = 1),
              non.zero.pro = checkNumeric(name.list$non.zero.pro, lower = 0, upper = 1),
              mean.coef0 = checkNumber(name.list$mean.coef0),
              mean.coef1 = checkNumber(name.list$mean.coef1),
              mean.coef2 = checkNumber(name.list$mean.coef2),
              noise.mean = checkNumber(name.list$noise.mean),
              noise.sd = checkNumber(name.list$noise.sd),
              sparse.fac = checkNumber(name.list$sparse.fac, lower = 0))

  if (all(checks == TRUE)) {
    valid <- TRUE
  } else {
    valid <- checks[checks != TRUE]
    valid <- paste(names(valid), valid, sep = ": ")
  }

  return(valid)
})
