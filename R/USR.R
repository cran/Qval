#' Caculate under-specifcation rate (USR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug A The Q-matrix that has being validated.
#'
#' @details
#' The USR is defned as:
#' \deqn{
#'  USR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} > q_{ik}^{s})}{I × K}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \code{k}th attribute of item \code{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{s}} denotes \code{k}th attribute of item \code{i} in the suggested Q-matrix(\code{Q.sug}),
#' and \eqn{I(\cdot)} is the indicator function.
#'
#' @return
#' A numeric (USR index).
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' example.Q1 <- sim.Q(5, 30)
#' example.Q2 <- sim.MQ(example.Q1, 0.1)
#' USR <- getUSR(example.Q1, example.Q2)
#' print(USR)
#'
#' @export
#'
getUSR <- function(Q.true, Q.sug) {
  USR <- 0
  for(i in 1:nrow(Q.true))
    for(j in 1:ncol(Q.true)) {
      if(Q.true[i, j] > Q.sug[i, j]){
        USR <- USR + 1
      }
    }
  USR <- USR / (nrow(Q.true)*ncol(Q.true))
  return(USR)
}
