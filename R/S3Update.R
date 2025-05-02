#' @title Update Method for Various Objects
#' @description
#' The \code{update} function provides a unified interface for refreshing or modifying
#' existing analysis objects produced by the Qval package, including \code{\link[Qval]{CDM}},
#' \code{\link[Qval]{validation}}, and \code{\link[Qval]{sim.data}} classes. By passing additional arguments,
#' users can rerun fitting or simulation routines without reconstructing the entire object
#' from scratch.
#'
#' @param x An object of class \code{\link[Qval]{CDM}}, \code{\link[Qval]{validation}}, 
#'          and \code{\link[Qval]{sim.data}}.
#' @param ... Additional arguments specific to the method being updated:
#'   \itemize{
#'     \item For \code{CDM}: \code{Y}, \code{Q}, \code{model}, \code{method},
#'       \code{mono.constraint}, \code{maxitr}, \code{verbose}.
#'     \item For \code{validation}: \code{Y}, \code{Q}, \code{CDM.obj}, \code{par.method},
#'       \code{mono.constraint}, \code{model}, \code{method}, \code{search.method},
#'       \code{iter.level}, \code{maxitr}, \code{eps}, \code{alpha.level}, \code{criter},
#'       \code{verbose}.
#'     \item For \code{sim.data}: \code{Q}, \code{N}, \code{IQ}, \code{model},
#'       \code{distribute}, \code{control}, \code{verbose}.
#'   }
#'
#' @return An updated object of the same class as \code{x}, reflecting any changes
#'   in the supplied parameters.
#'
#' @details
#' The \code{update} methods internally extract the original call arguments
#' from the input object, combine them with any new parameters provided in
#' \code{...}, and re-invoke the corresponding constructor (\code{\link[Qval]{CDM}}, 
#' \code{\link[Qval]{validation}}, and \code{\link[Qval]{sim.data}}). This approach ensures consistency
#' and preserves all untouched settings from the original object.
#'
#' @examples
#' set.seed(123)
#' library(Qval)
#' 
#' \donttest{
#' ################################################################
#' # Example 1: summary a CDM object                              #
#' ################################################################
#' Q <- sim.Q(3, 20)
#' IQ <- list(P0 = runif(20, 0, 0.2), P1 = runif(20, 0.8, 1))
#' data.obj <- sim.data(Q, N = 500, IQ = IQ, 
#'                      model = "GDINA", distribute = "horder")
#' CDM.obj <- CDM(data.obj$dat, Q, model = "GDINA", method = "BM")
#' summary(CDM.obj)
#' 
#' CDM.updated <- update(CDM.obj, method = "EM", maxitr = 1000)
#' summary(CDM.updated)
#' 
#'
#' ################################################################
#' # Example 2: summary a validation object                       #
#' ################################################################
#' MQ <- sim.MQ(Q, 0.1)
#' CDM.obj2 <- CDM(data.obj$dat, MQ)
#' validation.obj <- validation(data.obj$dat, MQ, CDM.obj2, 
#'                              method = "GDI")
#' summary(validation.obj)
#' 
#' validation.updated <- update(validation.obj, method = "Hull")
#' summary(validation.updated)
#' 
#'
#' ################################################################
#' # Example 3: summary a sim.data object                         #
#' ################################################################
#' data.obj2 <- sim.data(Q = sim.Q(3, 10), N = 1000)
#' summary(data.obj2)
#' 
#' data.updated <- update(data.obj2, N = 200)
#' summary(data.updated)
#' }
#'
#' @name update
update <- function(x, ...) {
  UseMethod("update")
}

#' @describeIn update Update method for CDM objects
#' @importFrom utils modifyList
#' @export
update.CDM <- function(x, ...) {
  arguments.current <- list(
    Y = x$arguments$Y,
    Q = x$arguments$Q,
    model = x$arguments$model,
    method = x$arguments$method,
    mono.constraint = x$arguments$mono.constraint,
    maxitr = x$arguments$maxitr,
    verbose = x$arguments$verbose
  )
  
  arguments.new <- modifyList(arguments.current, list(...))
  
  CDM.updated <- CDM(
    Y = arguments.new$Y,
    Q = arguments.new$Q,
    model = arguments.new$model,
    method = arguments.new$method,
    mono.constraint = arguments.new$mono.constraint,
    maxitr = arguments.new$maxitr,
    verbose = arguments.new$verbose
  )
  return(CDM.updated)
}

#' @describeIn update Update method for validation objects
#' @importFrom utils modifyList
#' @export
update.validation <- function(x, ...) {
  arguments.current <- list(
    Y = x$arguments$Y,
    Q = x$arguments$Q,
    CDM.obj = x$arguments$CDM.obj,
    par.method = x$arguments$par.method,
    mono.constraint = x$arguments$mono.constraint,
    model = x$arguments$model,
    method = x$arguments$method,
    search.method = x$arguments$search.method,
    iter.level = x$arguments$iter.level,
    maxitr = x$arguments$maxitr,
    eps = x$arguments$eps,
    alpha.level = x$arguments$alpha.level,
    criter = x$arguments$criter,
    verbose = x$arguments$verbose
  )
  
  arguments.new <- modifyList(arguments.current, list(...))
  
  validation.updated <- validation(
    Y = arguments.new$Y,
    Q = arguments.new$Q,
    CDM.obj = arguments.new$CDM.obj,
    par.method = arguments.new$par.method,
    mono.constraint = arguments.new$mono.constraint,
    model = arguments.new$model,
    method = arguments.new$method,
    search.method = arguments.new$search.method,
    iter.level = arguments.new$iter.level,
    maxitr = arguments.new$maxitr,
    eps = arguments.new$eps,
    alpha.level = arguments.new$alpha.level,
    criter = arguments.new$criter,
    verbose = arguments.new$verbose
  )
  return(validation.updated)
}

#' @describeIn update Update method for sim.data objects
#' @importFrom utils modifyList
#' @export
update.sim.data <- function(x, ...) {
  arguments.current <- list(
    Q = x$arguments$Q,
    N = x$arguments$N,
    IQ = x$arguments$IQ,
    model = x$arguments$model,
    distribute = x$arguments$distribute,
    control = x$arguments$control,
    verbose = x$arguments$verbose
  )
  arguments.new <- modifyList(arguments.current, list(...))
  
  data.updated <- sim.data(
    Q = arguments.new$Q,
    N = arguments.new$N,
    IQ = arguments.new$IQ,
    model = arguments.new$model,
    distribute = arguments.new$distribute,
    control = arguments.new$control,
    verbose = arguments.new$verbose
  )
  return(data.updated)
}
