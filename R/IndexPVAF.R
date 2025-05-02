#' Calculate PVAF
#'
#' @description
#' The function is able to calculate the proportion of variance accounted for (\eqn{PVAF}) for all items
#' after fitting \code{\link[Qval]{CDM}} or directly.
#'
#' @param Y A required \eqn{N} × \eqn{I} matrix or \code{data.frame} consisting of the responses of \code{N} individuals
#'          to \eqn{N} × \eqn{I} items. Missing values need to be coded as \code{NA}.
#' @param Q A required binary \eqn{I} × \eqn{K} matrix containing the attributes not required or required 
#'          master the items. The \code{i}th row of the matrix is a binary indicator vector indicating which
#'          attributes are not required (coded by 0) and which attributes are required (coded by 1) to master
#'          item \eqn{i}.
#' @param CDM.obj An object of class \code{CDM.obj}. Can be \code{NULL}, but when it is not \code{NULL}, it enables
#'                rapid verification of the Q-matrix without the need for parameter estimation.
#' @param model Type of model to be fitted; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"},
#'              \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#'
#' @details
#'  The intrinsic essence of the GDI index (as denoted by \eqn{\zeta_{2}}) is the weighted variance of
#'  all \eqn{2^{K\ast}} attribute mastery patterns' probabilities of correctly responding to
#'  item \eqn{i}, which can be computed as:
#'  \deqn{
#'  \zeta^2 =
#'  \sum_{l=1}^{2^K} \pi_{l}{(P(X_{pi}=1|\boldsymbol{\alpha}_{l}) - \bar{P}_{i})}^2
#' }
#' where \eqn{\pi_{l}} represents the prior probability of mastery pattern \eqn{l};
#' \eqn{\bar{P}_{i}=\sum_{k=1}^{2^K}\pi_{l}P(X_{pi}=1|\boldsymbol{\alpha}_{l})} is the weighted average of the correct
#' response probabilities across all attribute mastery patterns. When the q-vector
#' is correctly specified, the calculated \eqn{\zeta^2} should be maximized, indicating
#' the maximum discrimination of the item.
#'
#' Theoretically, \eqn{\zeta^{2}} is larger when \eqn{\boldsymbol{q}_{i}} is either specified correctly or over-specified,
#' unlike when \eqn{\boldsymbol{q}_{i}} is under-specified, and that when \eqn{\boldsymbol{q}_{i}} is over-specified, \eqn{\zeta^{2}}
#' is larger than but close to the value of \eqn{\boldsymbol{q}_{i}} when specified correctly. The value of \eqn{\zeta^{2}} continues to
#' increase slightly as the number of over-specified attributes increases, until \eqn{\boldsymbol{q}_{i}} becomes \eqn{\boldsymbol{q}_{i1:K}} 
#' (\eqn{\boldsymbol{q}_{i1:K}} = [11...1]).
#' Thus, \eqn{\zeta^{2} / \zeta_{max}^{2}} is computed to indicate the proportion of variance accounted for by \eqn{\boldsymbol{q}_{i}}, 
#' called the \eqn{PVAF}.
#'
#' @seealso \code{\link[Qval]{validation}}
#'
#' @return
#' An object of class \code{matrix}, which consisted of \eqn{PVAF} for each item and each possible q-vector.
#'
#' @author
#' Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' de la Torre, J., & Chiu, C. Y. (2016). A General Method of Empirical Q-matrix Validation. Psychometrika, 81(2), 253-273. DOI: 10.1007/s11336-015-9467-8.
#'
#' @examples
#' \donttest{
#' library(Qval)
#'
#' set.seed(123)
#'
#' ## generate Q-matrix and data
#' K <- 3
#' I <- 20
#' Q <- sim.Q(K, I)
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' data <- sim.data(Q = Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")
#'
#' ## calculate PVAF directly
#' PVAF <-get.PVAF(Y = data$dat, Q = Q)
#' print(PVAF)
#'
#' ## calculate PVAF after fitting CDM
#' CDM.obj <- CDM(data$dat, Q, model="GDINA")
#' PVAF <-get.PVAF(CDM.obj = CDM.obj)
#' print(PVAF)
#' }
#'
#' @export
#' @importFrom GDINA attributepattern
#'

get.PVAF <- function(Y = NULL, Q = NULL, CDM.obj = NULL, model = "GDINA"){

  if(is.null(CDM.obj) & (is.null(Y) | is.null(Q)))
    stop("one of [CDM.obj)] and [Y and Q] must not be NULL !!!")

  if(!is.null(Q)){
    I <- nrow(Q)
    K <- ncol(Q)
    L <- 2^K
    N <- nrow(Y)
  }else{
    I <- length(CDM.obj$analysis.obj$catprob.parm)
    K <- log2(length(CDM.obj$P.alpha))
    L <- length(CDM.obj$P.alpha)
    Y <- CDM.obj$analysis.obj$Y
    Q <- CDM.obj$analysis.obj$Q
    N <- nrow(Y)
  }

  # Create PVAF matrix
  PVAF <- matrix(NA, I, L)
  
  # Generate pattern names efficiently using apply and paste
  pattern <- attributepattern(K)
  pattern.names <- apply(pattern, 1, paste0, collapse = "")
  
  # Assign column and row names
  colnames(PVAF) <- pattern.names
  if(!is.null(rownames(Q))){
    rownames(PVAF) <- rownames(Q)
  }else{
    rownames(PVAF) <- paste0("item ", 1:I)
  }

  if(is.null(CDM.obj))
      CDM.obj <- CDM(Y=Y, Q=Q, model=model)
  alpha.P <- CDM.obj$alpha.P
  P.alpha <- CDM.obj$P.alpha
  alpha <- CDM.obj$alpha
  P.alpha.Xi <- CDM.obj$P.alpha.Xi

  for(i in 1:I){
    P.est <- calculatePEst(Y[, i], P.alpha.Xi)
    P.mean <- sum(P.est * P.alpha)

    P.Xi.alpha <- sapply(2:L, function(l) P_GDINA(pattern[l, ], P.est, pattern, P.alpha))
    zeta2 <- c(0, apply(P.Xi.alpha, 2, function(x) sum((x - P.mean)^2 * P.alpha)))
    PVAF[i, ] <- zeta2 / zeta2[L]
  }

  return(PVAF)
}
