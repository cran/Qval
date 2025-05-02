
#' 
#' Wald Test for Two Q-vectors
#' 
#' @description
#' This function flexibly provides the Wald test for any two q-vectors of a given item in the Q-matrix, 
#' but requires that the two q-vectors differ by only one attribute. Additionally, this function needs 
#' to accept a \code{CDM.obj}.
#' 
#' @details
#' \deqn{
#'    Wald = \left[\boldsymbol{R} \times \boldsymbol{P}_{i}(\boldsymbol{\alpha})\right]^{'}
#'    (\boldsymbol{R} \times \boldsymbol{V}_{i} \times \boldsymbol{R})^{-1}
#'    \left[\boldsymbol{R} \times P_{i}(\boldsymbol{\alpha})\right]
#' }
#' 
#' @param CDM.obj An object of class \code{CDM.obj}. @seealso \code{\link[Qval]{CDM}}.
#' @param q1 A q-vector
#' @param q2 Another q-vector
#' @param i the item needed to be validated
#'
#' @returns An object of class \code{htest} containing the following components:
#' \describe{
#'  \item{statistic}{The statistic of the Wald test.}
#'  \item{parameter}{the degrees of freedom for the Wald-statistic.}
#'  \item{p.value}{The p value}
#' }
#'
#' @examples
#' library(Qval)
#' set.seed(123)
#' 
#' K <- 3
#' I <- 20
#' N <- 500
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' Q <- sim.Q(K, I)
#' data <- sim.data(Q = Q, N = N, IQ = IQ, model = "GDINA", distribute = "horder")
#' 
#' CDM.obj <- CDM(data$dat, Q)
#' 
#' q1 <- c(1, 0, 0)
#' q2 <- c(1, 1, 0)
#' 
#' ## Discuss whether there is a significant difference when 
#' ## the q-vector of the 2nd item in the Q-matrix is q1 or q2.
#' Wald.test.obj <- Wald.test(CDM.obj, q1, q2, i=2)
#' 
#' print(Wald.test.obj)
#' 
#' 
#' @export
#' @import GDINA
#' @importFrom GDINA attributepattern LC2LG extract indlogPost score
#' @importFrom stats aggregate pchisq
#' @importFrom MASS ginv
#' 
Wald.test <- function(CDM.obj, q1, q2, i=1){
  
  if(sum(q1) > sum(q2)){
    temp <- q1
    q1 <- q2
    q2 <- temp
  }
  
  GDINA.obj <- CDM.obj$analysis.obj
  Y <- GDINA.obj$Y
  Q <- CDM.obj$analysis.obj$Q
  names.items <- rownames(Q)
  att.posi.i.k <- which(q2 > 0)
  att.posi.i <- which(q1 > 0)
  att.dif <- which(q1 != q2)
  
  Qr <- GDINA.obj$Q
  Qr[i, seq_len(ncol(Qr))] <- 0
  Qr[i, c(att.posi.i, att.dif)] <- 1
  etas <- LC2LG(as.matrix(Qr))
  expectedR <- GDINA::extract(GDINA.obj,"expectedCorrect.LC")
  expectedN <- GDINA::extract(GDINA.obj,"expectedTotal.LC")
  
  itemparj <- GDINA.obj$catprob.parm
  itemparj[[i]] <- aggregate(expectedR[i, ], by=list(etas[i, ]), sum)$x / aggregate(expectedN[i, ], by=list(etas[i, ]), sum)$x
  index <- data.frame(Cat=rep(1:length(rowSums(Qr) ),2^rowSums(Qr)) )
  index$Column <- seq_len(length(index$Cat))
  
  sco <- score(GDINA.obj, parm="prob") # a list with # of category elements
  sco[[i]] <- score_pj(Xj = Y[, i],                   # a vector of item responses to item i
                       parloc.j=etas[i, ,drop=FALSE],         # parameter locations for item i - H by 2^K matrix
                       catprob.j=itemparj[i],        # a list with H elements giving the reduced catprob.parm for each nonzero category
                       logpost=indlogPost(GDINA.obj))[[1]]
  
  ## v is the inversed information matrix
  if(GDINA::extract(GDINA.obj, "att.dist") != "saturated"){
    v <- inverse_crossprod(do.call(cbind, sco))
  }else{
    v <- inverse_crossprod(do.call(cbind,sco[-length(sco)]))
  }
  
  ## extrac parameters for item i
  cov <- v[index$Column[which(index$Cat == i)], index$Column[which(index$Cat == i)]]
  P.Xi.alpha.reduced <- itemparj[[i]]
  
  ## Restricted matrix based on q1 and q2
  Rmatrix <- get.Rmatrix(q1, q2)
  
  R.cov.Rt <- Rmatrix %*% cov %*% t(Rmatrix)  # R %*% cov %*% R^T
  R.P <- Rmatrix %*% P.Xi.alpha.reduced    # R %*% P.Xi.alpha.reduced
  
  R.cov.Rt.det <- det(R.cov.Rt)
  if (!is.na(R.cov.Rt.det)) {
    R.cov.Rt.solved <- solve(R.cov.Rt)
  }else{
    R.cov.Rt.solved <- ginv(R.cov.Rt)
  }
  Wald.statistic <- t(R.P) %*% R.cov.Rt.solved %*%R.P
  
  parameter = nrow(Rmatrix)
  p.value <- pchisq(Wald.statistic, parameter, lower.tail = FALSE)
  
  Wald.obj <- list(
    statistic = c(Wald = as.numeric(Wald.statistic)),
    parameter = c(df = parameter),
    p.value = as.numeric(p.value),
    method = "Wald test for two q-vectors",
    data.name = paste0(paste0("[", paste(q1, collapse = ""), "] "), 
                       "vs. ", 
                       paste0("[", paste(q2, collapse = ""), "] "), 
                       "for ", names.items[i])
  )
  class(Wald.obj) <- "htest"

  return(Wald.obj)
}
