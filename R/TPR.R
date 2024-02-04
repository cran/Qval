#' Caculate true-positive rate (TPR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.orig The Q-matrix need to be validated.
#' @param Q.sug The Q-matrix that has being validated.
#'
#' @details
#' TPR is defned as the proportion of correct elements which are correctly retained:
#' \deqn{
#'  TPR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} = q_{ik}^{s} | q_{ik}^{t} = q_{ik}^{o})}
#'  {\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} = q_{ik}^{o})}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \code{k}th attribute of item \eqn{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{o}} denotes \code{k}th attribute of item \code{i} in the original Q-matrix(\code{Q.orig}),
#' \eqn{q_{ik}^{s}} denotes \code{k}th attribute of item \code{i} in the suggested Q-matrix(\code{Q.sug}),
#' and \eqn{I(\cdot)} is the indicator function.
#'
#' @return
#' A numeric (TPR index).
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' example.Q1 <- sim.Q(5, 30)
#' example.Q2 <- sim.MQ(example.Q1, 0.1)
#' example.Q3 <- sim.MQ(example.Q1, 0.05)
#' TPR <- getTPR(example.Q1, example.Q2, example.Q3)
#'
#' print(TPR)
#'
#' @export
#'
getTPR <- function(Q.true, Q.orig, Q.sug) {
  TPR <- 0
  sum <- 0
  for(i in 1:nrow(Q.true))
    for(j in 1:ncol(Q.true)) {
      if(Q.true[i, j] == Q.orig[i, j]){
        sum <- sum + 1
        if(Q.sug[i, j] == Q.true[i, j])
          TPR <- TPR + 1
      }
    }
  TPR <- TPR / sum
  return(TPR)
}
