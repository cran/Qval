
#' 
#' Wald.test for two q-vecotrs
#' 
#' @description
#' This function flexibly provides the Wald test for any two q-vectors of a given item in the Q-matrix, 
#' but requires that the two q-vectors differ by only one attribute. Additionally, this function needs 
#' to accept a \code{CDM.obj}.
#' 
#' @details
#' \deqn{
#'    Wald = (\mathbf{R} \times P_{i}(\mathbf{\alpha}))^{'}
#'    (\mathbf{R} \times \mathbf{V}_{i} \times \mathbf{R})^{-1}
#'    (\mathbf{R} \times P_{i}(\mathbf{\alpha}))
#' }
#' 
#' @param CDM.obj An object of class \code{CDM.obj}. @seealso \code{\link[Qval]{CDM}}.
#' @param Q.i A q-vector
#' @param Q.i.k Another q-vector
#' @param i the item you focusing on
#'
#' @returns An object of class \code{list} containing the following components:
#' \item{Wald.statistic}{The statistic of the Wald test.}
#' \item{p.value}{The p value}
#'
#' @examples
#' 
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
#' Q.i <- c(1, 0, 0)
#' Q.i.k <- c(1, 1, 0)
#' 
#' ## Discuss whether there is a significant difference when 
#' ## the q-vector of the 2nd item in the Q-matrix is Q.i or Q.i.k.
#' Wald.test.obj <- Wald.test(CDM.obj, Q.i, Q.i.k, i=2)
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
Wald.test <- function(CDM.obj, Q.i, Q.i.k, i=1){
  
  GDINA.obj <- CDM.obj$analysis.obj
  Y <- GDINA.obj$Y
  att.posi.i.k <- which(Q.i.k > 0)
  att.posi.i <- which(Q.i > 0)
  att.dif <- which(Q.i != Q.i.k)
  
  Qr <- extract(GDINA.obj,"Q")
  Qr[i, seq_len(ncol(Qr))] <- 0
  Qr[i, c(att.posi.i, att.dif)] <- 1
  etas <- LC2LG(as.matrix(Qr))
  itemparj <- extract(GDINA.obj,"catprob.parm")
  expectedR <- extract(GDINA.obj,"expectedCorrect.LC")
  expectedN <- extract(GDINA.obj,"expectedTotal.LC")
  
  itemparj[[i]] <- aggregate(expectedR[i, ], by=list(etas[i, ]), sum)$x / aggregate(expectedN[i, ], by=list(etas[i, ]), sum)$x
  index <- data.frame(Cat=rep(1:length(rowSums(Qr) ),2^rowSums(Qr)) )
  index$Column <- seq_len(length(index$Cat))
  
  sco <- score(GDINA.obj, parm="prob") # a list with # of category elements
  sco[[i]] <- score_pj(Xj = Y[, i],                   # a vector of item responses to item i
                               parloc.j=etas[i, ,drop=FALSE],         # parameter locations for item i - H by 2^K matrix
                               catprob.j=itemparj[i],        # a list with H elements giving the reduced catprob.parm for each nonzero category
                               logpost=indlogPost(GDINA.obj))[[1]]
  
  ## v is the inversed information matrix
  if(extract(GDINA.obj, "att.dist") != "saturated"){
    v <- inverse_crossprod(do.call(cbind, sco))
  }else{
    v <- inverse_crossprod(do.call(cbind,sco[-length(sco)]))
  }
  
  ## extrac parameters for item i
  cov.cur <- v[index$Column[which(index$Cat == i)], index$Column[which(index$Cat == i)]]
  P.Xi.alpha.reduced <- itemparj[[i]]
  
  ## Restricted matrix based on Q.i and Q.i.k
  Rmatrix <- get.Rmatrix(Q.i, Q.i.k)
  
  Wald.statistic <-t(Rmatrix %*% P.Xi.alpha.reduced) %*% ginv(Rmatrix %*% cov.cur %*% t(Rmatrix)) %*% (Rmatrix %*% P.Xi.alpha.reduced)
  p.value <- pchisq(Wald.statistic, nrow(Rmatrix), lower.tail = FALSE)
  
  return(list(Wald.statistic=Wald.statistic, p.value=p.value))
}
