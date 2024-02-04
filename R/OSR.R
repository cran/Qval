#' Caculate over-specifcation rate (OSR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug The Q-matrix that has being validated.
#'
#' @details
#'
#' The OSR is defned as:
#' \deqn{
#'  OSR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} < q_{ik}^{s})}{I × K}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \code{k}th attribute of item \code{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{s}} denotes \code{k}th attribute of item \code{i} in the suggested Q-matrix(\code{Q.sug}),
#' and \eqn{I(\cdot)} is the indicator function.
#'
#' @return
#' A numeric (OSR index).
#'
#' @examples
#'
#' library(Qval)
#'
#' set.seed(123)
#'
#' example.Q1 <- sim.Q(5, 30)
#' example.Q2 <- sim.MQ(example.Q1, 0.1)
#' OSR <- getOSR(example.Q1, example.Q2)
#' print(OSR)
#'
#' @export
#'
getOSR <- function(Q.true, Q.sug) {
  OSR <- 0
  for(i in 1:nrow(Q.true))
    for(j in 1:ncol(Q.true)) {
      if(Q.true[i, j] < Q.sug[i, j]){
        OSR <- OSR + 1
      }
    }
  OSR <- OSR / (nrow(Q.true)*ncol(Q.true))
  return(OSR)
}
