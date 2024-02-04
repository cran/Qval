keep.convex <- function(Hull, fjk, npk){

  Hull.ori <- Hull

  K <- length(Hull) + 1
  Hull[which(Hull < 1)] <- 0
  while(any(Hull == 0) & all(!is.nan(Hull))){
    pre.posi <- nex.posi <- rep(0, K - 1)

    for(k in 1:(K-1)){
      if(all(!is.nan(Hull))){
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
          Hull[k] <- (pre+1e-10) / (nex+2e-10)
        }
        Hull[which(Hull == 0)] <- -Inf
        Hull[which(Hull < 1 & Hull != -Inf)] <- 0
      }
    }
  }

  if(any(is.nan(Hull)))
    Hull <- Hull.ori

  return(Hull)
}
