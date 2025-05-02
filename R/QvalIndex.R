#' Calculate Over-Specification Rate (OSR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug The Q-matrix that has been validated.
#'
#' @details
#'
#' The OSR is defned as:
#' \deqn{
#'  OSR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} < q_{ik}^{s})}{I × K}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \eqn{k}th attribute of item \eqn{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{s}} denotes \eqn{k}th attribute of item \eqn{i} in the suggested Q-matrix(\code{Q.sug}),
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
#' Q1 <- sim.Q(5, 30)
#' Q2 <- sim.MQ(Q1, 0.1)
#' OSR <- zOSR(Q1, Q2)
#' print(OSR)
#'
#' @export
#'
zOSR <- function(Q.true, Q.sug) {
  # Compare matrices directly element-wise and sum the logical values (TRUE = 1, FALSE = 0)
  OSR <- sum(Q.true < Q.sug) / (length(Q.true))  # length(Q.true) is equivalent to nrow(Q.true) * ncol(Q.true)
  return(OSR)
}


#' Calculate Q-matrix Recovery Rate (QRR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug The Q-matrix that has been validated.
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
#' Q1 <- sim.Q(5, 30)
#' Q2 <- sim.MQ(Q1, 0.1)
#' QRR <- zQRR(Q1, Q2)
#' print(QRR)
#'
#' @export
#'
zQRR <- function(Q.true, Q.sug) {
  return(1 - sum(abs(Q.true - Q.sug)) / prod(dim(Q.true)))
}


#' Calculate True-Negative Rate (TNR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.orig The Q-matrix need to be validated.
#' @param Q.sug The Q-matrix that has been validated.
#'
#' @details
#' TNR is defined as the proportion of correct elements which are correctly retained:
#' \deqn{
#'  TNR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} = q_{ik}^{s} | q_{ik}^{t} \neq q_{ik}^{o})}
#'  {\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} \neq q_{ik}^{o})}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \eqn{k}th attribute of item \eqn{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{o}} denotes \eqn{k}th attribute of item \eqn{i} in the original Q-matrix(\code{Q.orig}),
#' \eqn{q_{ik}^{s}} denotes \eqn{k}th attribute of item \eqn{i} in the suggested Q-matrix(\code{Q.sug}),
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
#' Q1 <- sim.Q(5, 30)
#' Q2 <- sim.MQ(Q1, 0.1)
#' Q3 <- sim.MQ(Q1, 0.05)
#' TNR <- zTNR(Q1, Q2, Q3)
#'
#' print(TNR)
#'
#' @export
#'
zTNR <- function(Q.true, Q.orig, Q.sug) {
  # Compare the true values with the original to find where they differ
  diff_indices <- Q.true != Q.orig
  
  # Calculate TNR where differences occur
  TNR <- sum(diff_indices & (Q.sug == Q.true)) / sum(diff_indices)
  
  return(TNR)
}


#' Calculate True-Positive Rate (TPR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.orig The Q-matrix need to be validated.
#' @param Q.sug The Q-matrix that has been validated.
#'
#' @details
#' TPR is defned as the proportion of correct elements which are correctly retained:
#' \deqn{
#'  TPR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} = q_{ik}^{s} | q_{ik}^{t} = q_{ik}^{o})}
#'  {\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} = q_{ik}^{o})}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \eqn{k}th attribute of item \eqn{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{o}} denotes \eqn{k}th attribute of item \eqn{i} in the original Q-matrix(\code{Q.orig}),
#' \eqn{q_{ik}^{s}} denotes \eqn{k}th attribute of item \eqn{i} in the suggested Q-matrix(\code{Q.sug}),
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
#' Q1 <- sim.Q(5, 30)
#' Q2 <- sim.MQ(Q1, 0.1)
#' Q3 <- sim.MQ(Q1, 0.05)
#' TPR <- zTPR(Q1, Q2, Q3)
#'
#' print(TPR)
#'
#' @export
#'
zTPR <- function(Q.true, Q.orig, Q.sug) {
  # Compare the true values with the original to find where they match
  match_indices <- Q.true == Q.orig
  
  # Calculate TPR where matches occur
  TPR <- sum(match_indices & (Q.sug == Q.true)) / sum(match_indices)
  
  return(TPR)
}


#' Calculate Under-Specification Rate (USR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug The Q-matrix that has been validated.
#'
#' @details
#' The USR is defned as:
#' \deqn{
#'  USR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} > q_{ik}^{s})}{I × K}
#' }
#' where \eqn{q_{ik}^{t}} denotes the \eqn{k}th attribute of item \eqn{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{q_{ik}^{s}} denotes \eqn{k}th attribute of item \eqn{i} in the suggested Q-matrix(\code{Q.sug}),
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
#' Q1 <- sim.Q(5, 30)
#' Q2 <- sim.MQ(Q1, 0.1)
#' USR <- zUSR(Q1, Q2)
#' print(USR)
#'
#' @export
#'
zUSR <- function(Q.true, Q.sug) {
  # Count how many times Q.true is greater than Q.sug
  USR <- sum(Q.true > Q.sug) / (nrow(Q.true) * ncol(Q.true))
  
  return(USR)
}


#' Calculate Vector Recovery Ratio (VRR)
#'
#' @param Q.true The true Q-matrix.
#' @param Q.sug The Q-matrix that has been validated.
#'
#' @details
#' The VRR shows the ability of the validation method to recover q-vectors, and is determined by
#' \deqn{
#'  VRR =\frac{\sum_{i=1}^{I}I(\mathbf{q}_{i}^{t} = \mathbf{q}_{i}^{s})}{I}
#' }
#' where \eqn{\mathbf{q}_{i}^{t}} denotes the \eqn{\mathbf{q}}-vector of item \eqn{i} in the true Q-matrix (\code{Q.true}),
#' \eqn{\mathbf{q}_{i}^{s}} denotes the \eqn{\mathbf{q}}-vector of item \eqn{i} in the suggested Q-matrix(\code{Q.sug}),
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
#' Q1 <- sim.Q(5, 30)
#' Q2 <- sim.MQ(Q1, 0.1)
#' VRR <- zVRR(Q1, Q2)
#' print(VRR)
#'
#' @export
#'
zVRR <- function(Q.true, Q.sug) {
  # Count how many rows are identical between Q.true and Q.sug
  same <- sum(rowSums(Q.true != Q.sug) == 0)
  
  # Calculate the ratio of identical rows
  return(same / nrow(Q.true))
}


