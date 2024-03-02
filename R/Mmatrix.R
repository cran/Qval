#' Calculate \eqn{\mathbf{M}} matrix
#'
#' @description
#' Calculate \eqn{\mathbf{M}} matrix for stauted CDMs (de la Torre, 2011).
#'
#' @param K The number of attributes. Can be NULL if \code{pattern} is passed to the function and is not NULL.
#' @param pattern The knowledge state matrix containing all possible attribute mastery pattern.
#'                Can be gained from @seealso \code{\link[GDINA]{attributepattern}}. Also can be NULL if \code{K}
#'                is passed to the function and is not NULL.
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
#' example.Mmatrix <-  get.Mmatrix(K = 5)
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
  Mmatrix <- matrix(0, L, L)
  for(l1 in 1:L)
    for(l2 in 1:L)
      Mmatrix[l1, l2] <- prod(pattern[l1, ]^pattern[l2, ])
  return(Mmatrix)
}
