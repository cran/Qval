#' Caculate Q-matrix recovery rate (QRR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug A The Q-matrix that has being validated.
#'
#' @details
#' The Q-matrix recovery rate (QRR) provides information on overall performance, and is defned as:
#' \deqn{
#'  QRR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} = q_{ik}^{s})}{I × K}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \eqn{k}th attribute of item \eqn{i} in the true Q-matrix (\eqn{Q.true}),
#' \eqn{q_{ik}^{s}} denotes \eqn{k}th attribute of item \eqn{i} in the suggested Q-matrix(\eqn{Q.sug}),
#' and \eqn{I(\cdot)} is the indicator function.
#'
#' @return
#' A numeric (QRR index).
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' example.Q1 <- sim.Q(5, 30)
#' example.Q2 <- sim.MQ(example.Q1, 0.1)
#' QRR <- getQRR(example.Q1, example.Q2)
#' print(QRR)
#'
#' @export
#'
getQRR <- function(Q.true, Q.sug) {
  return(1 - sum(abs(Q.true - Q.sug)) / (nrow(Q.true)*ncol(Q.true)))
}
