#' Calculate \eqn{\mathbf{M}} matrix
#'
#' @description
#' Calculate \eqn{\mathbf{M}} matrix for saturated CDMs (de la Torre, 2011). The \eqn{\mathbf{M}} matrix is a matrix used to 
#' represent the interaction mechanisms between attributes.  
#'
#' @param K The number of attributes. Can be \code{NULL} if the argument \code{pattern} is not \code{NULL}.
#' @param pattern The attribute mastery pattern matrix containing all possible attribute mastery pattern.
#'                Can be gained from \code{\link[GDINA]{attributepattern}}. Also can be \code{NULL} if \code{K}
#'                is not \code{NULL}.
#'
#' @return
#' An object of class \code{matrix}.
#'
#' @author
#' Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' de la Torre, J. (2011). The Generalized DINA Model Framework. Psychometrika, 76(2), 179-199. DOI: 10.1007/s11336-011-9207-7.
#'
#' @examples
#'
#' library(Qval)
#'
#' Mmatrix <-  get.Mmatrix(K = 3)
#' 
#' print(Mmatrix)
#'
#' @export
#' @importFrom GDINA attributepattern
#'

get.Mmatrix <- function(K = NULL, pattern = NULL){

  if(is.null(K) & is.null(pattern)){
    stop("one of K and pattern must not be NULL !!!")
  }else if(is.null(pattern)){
    pattern <- attributepattern(K)
  }

  L <- nrow(pattern)
  K <- ncol(pattern)

  Mmatrix <- outer(1:L, 1:L, Vectorize(function(l1, l2) prod(pattern[l1, ]^pattern[l2, ])))
  
  return(Mmatrix)
}
