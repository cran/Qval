
#' 
#' Restriction Matrix
#' 
#' @description
#' This function returns the restriction matrix (de la Torre, 2011; Ma & de la Torre, 2020) based on two q-vectors, 
#' where the two q-vectors can only differ by one attribute.
#' 
#' @param q1 A q-vector
#' @param q2 Another q-vector
#' 
#' @seealso \code{\link[Qval]{Wald.test}}
#' 
#' @return A restriction matrix
#' 
#' @references
#' de la Torre, J. (2011). The Generalized DINA Model Framework. Psychometrika, 76(2), 179-199. DOI: 10.1007/s11336-011-9207-7.
#' 
#' Ma, W., & de la Torre, J. (2020). An empirical Q-matrix validation method for the sequential generalized DINA model. British Journal of Mathematical and Statistical Psychology, 73(1), 142-163. DOI: 10.1111/bmsp.12156.
#' 
#' @examples
#' q1 <- c(1, 1, 0)
#' q2 <- c(1, 1, 1)
#' 
#' Rmatrix <- get.Rmatrix(q1, q2)
#' 
#' print(Rmatrix)
#' 
#' 
#' @export
#' @importFrom GDINA attributepattern
#' 
get.Rmatrix <- function(q1, q2){
  att.posi.i.k <- which(q2 > 0)
  att.posi.i <- which(q1 > 0)
  
  pattern.reduced.q1 <- attributepattern(length(att.posi.i))  ## reduced KS pattern for q2
  pattern.reduced.q2 <- attributepattern(length(att.posi.i.k))  ## reduced KS pattern for q2
  
  Rmatrix <- matrix(0, nrow(pattern.reduced.q1), nrow(pattern.reduced.q2))   ## the Restricted matrix
  dif.reduced <- which(q1[att.posi.i.k] != q2[att.posi.i.k])  ## find the different position between q1 and q2
  pattern.compare <- matrix(pattern.reduced.q2[, -dif.reduced], nrow = nrow(pattern.reduced.q2))  ## Used to determine which two vectors need to be compared
  is.compared <- rep(0, nrow(pattern.compare))  ## Used to record which vectors have already been recored in Rmatrix
  
  t <- 1
  while (any(!is.compared)) {
    
    ## Start searching from vectors that have never been recorded
    posi <- which(is.compared == 0)[1]  
    
    for(l in c(1:nrow(pattern.compare))[-posi]){
      
      ## If true, it indicates that the vectors at positions posi and l need to be recorded in the Rmatrix.
      if(all(pattern.compare[posi, ] == pattern.compare[l, ])){
        is.compared[c(posi, l)] <- 1
        Rmatrix[t, sort(c(posi, l))] <- c(1, -1)
        t <- t + 1
      }
      
    }
  }
  
  return(Rmatrix)
}