#' @title Summary Methods for Various Objects
#'
#' @description
#' Generate concise summary statistics for objects created by the Qval package.
#' The output is a named list tailored to the class of the input:
#' \describe{
#'   \item{\code{\link[Qval]{CDM}}}{contains the original call, dataset dimensions, model fit object, and attribute-pattern distribution.}
#'   \item{\code{\link[Qval]{validation}}}{contains the original call, suggested Q-matrix, and original Q-matrix.}
#'   \item{\code{\link[Qval]{sim.data}}}{contains the original call, dataset dimensions, and attribute-pattern distribution.}
#' }
#'
#' @details
#' \describe{
#'   \item{call}{A string capturing the original function invocation.}
#'   \item{base}{A numeric vector \code{c(N, I, K)} giving the number of examinees (\eqn{N}), 
#'               the number of items (\eqn{I}), and the number of attributes (\eqn{K}).}
#'   \item{model.fit}{(CDM only) The fitted model object returned by \code{\link[Qval]{CDM}}.}
#'   \item{patterns}{(CDM and sim.data) A data.frame of frequencies (\code{freq}) and proportions 
#'                   (\code{prop}) of each attribute pattern.}
#'   \item{Q.sug}{(validation only) Suggested Q-matrix from \code{\link[Qval]{validation}}.}
#'   \item{Q.orig}{(validation only) Original Q-matrix provided to \code{\link[Qval]{sim.data}}.}
#' }
#'
#' @param object An object of class \code{\link[Qval]{CDM}}, \code{\link[Qval]{validation}}, or \code{\link[Qval]{sim.data}}.
#' @param ...   Currently unused. Additional arguments are ignored.
#'
#' @return A named list with class \code{summary.<class>} containing the components above.
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
#' CDM.obj <- CDM(data.obj$dat, Q, model = "GDINA", method = "EM")
#' summary(CDM.obj)
#' 
#'
#' ################################################################
#' # Example 2: summary a validation object                       #
#' ################################################################
#' MQ <- sim.MQ(Q, 0.1)
#' CDM.obj2 <- CDM(data.obj$dat, MQ)
#' val.obj <- validation(data.obj$dat, MQ, CDM.obj2, method = "GDI")
#' summary(val.obj)
#' 
#'
#' ################################################################
#' # Example 3: summary a sim.data object                         #
#' ################################################################
#' data.obj2 <- sim.data(Q = sim.Q(3, 10), N = 1000)
#' summary(data.obj2)
#' }
#' 
#' @name summary
NULL

#' @describeIn summary Summary method for CDM objects
#' @export
summary.CDM <- function(object, ...) {
  dat <- object$analysis.obj$Y
  Q <- object$analysis.obj$Q
  catprob.parm <- object$analysis.obj$catprob.parm
  alpha <- object$alpha
  model.fit <- object$model.fit
  
  call <- paste(deparse(extract.CDM(object,"call"), width.cutoff = 30), sep = "\n", collapse = "\n")
  
  base <- c(N=nrow(dat), I=ncol(dat), K=ncol(Q))
  
  alpha <- apply(alpha, 1, function(x) paste(x, collapse=""))
  alpha.freq <- table(alpha)
  alpha.prop <- round(alpha.freq / sum(alpha.freq), 4)
  patterns <- data.frame(rbind(as.character(alpha.freq), as.character(alpha.prop)))
  colnames(patterns) <- names(alpha.freq)
  rownames(patterns) <- c("freq", "prop")
  
  out <- list(
    call = call,
    base=base,
    model.fit = model.fit,
    patterns = patterns
  )
  
  class(out) <- "summary.CDM"
  return(out)
}

#' @describeIn summary Summary method for validation objects
#' @export
summary.validation <- function(object, ...){
  
  call <- paste(deparse(extract.CDM(object,"call"), width.cutoff = 30), sep = "\n", collapse = "\n")
  
  Q.sug <- data.frame(extract.validation(object, "Q.sug"))
  Q.orig <- data.frame(extract.validation(object, "Q.orig"))
  out <- list(
    call = call,
    Q.sug=Q.sug,
    Q.orig = Q.orig
  )
  
  class(out) <- "summary.validation"
  return(out)
}

#' @describeIn summary Summary method for sim.data objects
#' @export
summary.sim.data <- function(object, ...) {
  
  dat <- object$dat
  Q <- object$Q
  catprob.parm <- object$catprob.parm
  alpha <- object$attribute
  
  call <- paste(deparse(extract.CDM(object,"call"), width.cutoff = 30), sep = "\n", collapse = "\n")
  
  base <- c(N=nrow(dat), I=ncol(dat), K=ncol(Q))
  
  alpha <- apply(alpha, 1, function(x) paste(x, collapse=""))
  alpha.freq <- table(alpha)
  alpha.prop <- round(alpha.freq / sum(alpha.freq), 5)
  patterns <- data.frame(rbind(as.character(alpha.freq), as.character(alpha.prop)))
  colnames(patterns) <- names(alpha.freq)
  rownames(patterns) <- c("freq", "prop")
  
  out <- list(
    call = call,
    base=base,
    patterns = patterns
  )
  
  class(out) <- "summary.sim.data"
  return(out)
}