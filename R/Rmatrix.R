
#' 
#' Restriction matrix
#' 
#' @description
#' This function returns the restriction matrix (de la Torre, 2011; Ma & de la Torre, 2020) based on two q-vectors, 
#' where the two q-vectors can only differ by one attribute.
#' 
#' @param Q.i A q-vector
#' @param Q.i.k Another q-vector
#' 
#' @return A restriction matrix
#' 
#' @references
#' de la Torre, J. (2011). The Generalized DINA Model Framework. Psychometrika, 76(2), 179-199. DOI: 10.1007/s11336-011-9207-7.
#' 
#' Ma, W., & de la Torre, J. (2020). An empirical Q-matrix validation method for the sequential generalized DINA model. British Journal of Mathematical and Statistical Psychology, 73(1), 142-163. DOI: 10.1111/bmsp.12156.
#' 
#' @examples
#' Q.i <- c(1, 1, 0)
#' Q.i.k <- c(1, 1, 1)
#' 
#' Rmatrix <- get.Rmatrix(Q.i, Q.i.k)
#' 
#' print(Rmatrix)
#' 
#' 
#' @export
#' @importFrom GDINA attributepattern
#' 
get.Rmatrix <- function(Q.i, Q.i.k){
  att.posi.i.k <- which(Q.i.k > 0)
  att.posi.i <- which(Q.i > 0)
  
  pattern.reduced.Q.i <- attributepattern(length(att.posi.i))  ## reduced KS pattern for Q.i.k
  pattern.reduced.Q.i.k <- attributepattern(length(att.posi.i.k))  ## reduced KS pattern for Q.i.k
  
  Rmatrix <- matrix(0, nrow(pattern.reduced.Q.i), nrow(pattern.reduced.Q.i.k))   ## the Restricted matrix
  dif.reduced <- which(Q.i[att.posi.i.k] != Q.i.k[att.posi.i.k])  ## find the different position between Q.i and Q.i.k
  pattern.compare <- matrix(pattern.reduced.Q.i.k[, -dif.reduced], nrow = nrow(pattern.reduced.Q.i.k))  ## Used to determine which two vectors need to be compared
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