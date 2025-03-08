
#' 
#' the Wald test for two q-vectors
#' 
#' @description
#' This function flexibly provides the Wald test for any two q-vectors of a given item in the Q-matrix, 
#' but requires that the two q-vectors differ by only one attribute. Additionally, this function needs 
#' to accept a \code{CDM.obj}.
#' 
#' @details
#' \deqn{
#'    Wald = \left[\mathbf{R} \times P_{i}(\mathbf{\alpha})\right]^{'}
#'    (\mathbf{R} \times \mathbf{V}_{i} \times \mathbf{R})^{-1}
#'    \left[\mathbf{R} \times P_{i}(\mathbf{\alpha})\right]
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
  att.posi.i.k <- which(q2 > 0)
  att.posi.i <- which(q1 > 0)
  att.dif <- which(q1 != q2)
  
  Qr <- GDINA.obj$Q
  Qr[i, seq_len(ncol(Qr))] <- 0
  Qr[i, c(att.posi.i, att.dif)] <- 1
  etas <- LC2LG(as.matrix(Qr))
  itemparj <- GDINA.obj$catprob.parm
  expectedR <- GDINA::extract(GDINA.obj,"expectedCorrect.LC")
  expectedN <- GDINA::extract(GDINA.obj,"expectedTotal.LC")
  
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
  
  Wald.statistic <-t(Rmatrix %*% P.Xi.alpha.reduced) %*% ginv(Rmatrix %*% cov %*% t(Rmatrix)) %*% (Rmatrix %*% P.Xi.alpha.reduced)
  parameter = nrow(Rmatrix)
  p.value <- pchisq(Wald.statistic, parameter, lower.tail = FALSE)
  
  Wald.obj <- list(statistic=Wald.statistic, parameter=parameter, p.value=p.value)
  
  class(Wald.obj) <- "htest"
  
  return(Wald.obj)
}
