

# This function is used to maintain the convexity of the Hull plot. 
# The principle is to iteratively remove all points where the ratio of slopes is less than 1.
keep.convex <- function(Hull, fjk, npk){

  Hull.ori <- Hull

  K <- length(Hull) + 1
  Hull[which(Hull < 1)] <- 0
  
  if(any(Hull == 0)){
    while(any(Hull == 0) & all(!is.nan(Hull))){
      pre.posi <- nex.posi <- rep(0, K - 1)
      
      posi <- list()
      
      for(k in 1:(K-1)){
        if(all(!is.nan(Hull))){
          
          ## Determine the position of the point preceding point k.
          if(Hull[k] >= 1){
            if(k == 1){
              pre.posi[k] <- 1
            }else{
              for(kk in (k-1):1){
                if(Hull[kk] >= 1){
                  pre.posi[k] <- kk + 1
                  break
                }
                if(Hull[kk] < 1 & kk == 1)
                  pre.posi[k] <- 1
              }
            }
            
            ## Determine the position of the point following point k.
            if(k == (K-1)){
              nex.posi[k] <- K + 1
            }else{
              for(kk in (k+1):(K-1)){
                if(Hull[kk] >= 1){
                  nex.posi[k] <- kk + 1
                  break
                }
                if(Hull[kk] < 1 & kk == (K-1))
                  nex.posi[k] <- K+1
              }
            }
            pre <- (fjk[k+1]-fjk[pre.posi[k]]) / (npk[k+1]- npk[pre.posi[k]])
            nex <- (fjk[nex.posi[k]]-fjk[k+1]) / (npk[nex.posi[k]]- npk[k+1])
            Hull[k] <- (pre+1e-10) / (nex+1e-10)
            
            posi[[k]] <- c(pre.posi[k], k+1, nex.posi[k])
          }
          Hull[which(Hull == 0)] <- -Inf
          Hull[which(Hull < 1 & Hull != -Inf)] <- 0
        }
      }
    }
  }else{
    posi <- list()
    posi[[1]] <- 1:length(fjk)
  }

  if(any(is.nan(Hull)))
    Hull <- Hull.ori
  if(length(posi) == 0)
    posi <- list(c(1, length(fjk)))

  return(list(Hull=Hull, posi=posi))
}
