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
#' OSR <- zOSR(example.Q1, example.Q2)
#' print(OSR)
#'
#' @export
#'
zOSR <- function(Q.true, Q.sug) {
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
#' QRR <- zQRR(example.Q1, example.Q2)
#' print(QRR)
#'
#' @export
#'
zQRR <- function(Q.true, Q.sug) {
  return(1 - sum(abs(Q.true - Q.sug)) / (nrow(Q.true)*ncol(Q.true)))
}


#' Calculate true negative rate (TNR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.orig The Q-matrix need to be validated.
#' @param Q.sug The Q-matrix that has being validated.
#'
#' @details
#' TNR is defined as the proportion of correct elements which are correctly retained:
#' \deqn{
#'  TNR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} = q_{ik}^{s} | q_{ik}^{t} \neq q_{ik}^{o})}
#'  {\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} \neq q_{ik}^{o})}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \code{k}th attribute of item \code{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{o}} denotes \code{k}th attribute of item \code{i} in the original Q-matrix(\code{Q.orig}),
#' \eqn{q_{ik}^{s}} denotes \code{k}th attribute of item \code{i} in the suggested Q-matrix(\code{Q.sug}),
#' and \eqn{I(\cdot)} is the indicator function.
#'
#' @return
#' A numeric (TNR index).
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' example.Q1 <- sim.Q(5, 30)
#' example.Q2 <- sim.MQ(example.Q1, 0.1)
#' example.Q3 <- sim.MQ(example.Q1, 0.05)
#' TNR <- zTNR(example.Q1, example.Q2, example.Q3)
#'
#' print(TNR)
#'
#' @export
#'
zTNR <- function(Q.true, Q.orig, Q.sug) {
  TNR <- 0
  sum <- 0
  for(i in 1:nrow(Q.true))
    for(j in 1:ncol(Q.true)) {
      if(Q.true[i, j] != Q.orig[i, j]){
        sum <- sum + 1
        if(Q.sug[i, j] == Q.true[i, j])
          TNR <- TNR + 1
      }
    }
  TNR <- TNR / sum
  return(TNR)
}


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
#' TPR <- zTPR(example.Q1, example.Q2, example.Q3)
#'
#' print(TPR)
#'
#' @export
#'
zTPR <- function(Q.true, Q.orig, Q.sug) {
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
#' USR <- zUSR(example.Q1, example.Q2)
#' print(USR)
#'
#' @export
#'
zUSR <- function(Q.true, Q.sug) {
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



#' Caculate vector recovery ratio (VRR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug A The Q-matrix that has being validated.
#'
#' @details
#' The VRR shows the ability of the validation method to recover q-vectors, and is determined by
#' \deqn{
#'  VRR =\frac{\sum_{i=1}^{I}I(\mathbf{q}_{i}^{t} = \mathbf{q}_{i}^{s})}{I}
#' }
#' where \eqn{\mathbf{q}_{i}^{t}} denotes the \eqn{\mathbf{q}}-vector of item \code{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{\mathbf{q}_{i}^{s}} denotes the \eqn{\mathbf{q}}-vector of item \code{i} in the suggested Q-matrix(\code{Q.sug}),
#' and \eqn{I(\cdot)} is the indicator function.
#'
#' @return
#' A numeric (VRR index).
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' example.Q1 <- sim.Q(5, 30)
#' example.Q2 <- sim.MQ(example.Q1, 0.1)
#' VRR <- zVRR(example.Q1, example.Q2)
#' print(VRR)
#'
#' @export
#'
zVRR <- function(Q.true, Q.sug) {
  VRR <- 0
  K <- ncol(Q.true)
  I <- nrow(Q.true)
  dif <- rowSums(abs(Q.true - Q.sug))
  same <- length(which(dif == 0))
  return(same/I)
}

