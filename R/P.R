P.GDINA <- function(Qi, P.est, pattern, P.alpha){

  K <- length(Qi)
  L <- 2^K
  P.Xj.alpha <- rep(0, L)

  att <- which(Qi == 1)
  pattern.temp <- matrix(pattern[, att], nrow = L)
  pattern.paste <- as.character(pattern.temp[, 1])
  K <- ncol(pattern.temp)
  if(K > 1)
    for(k in 2:K)
      pattern.paste <- paste0(pattern.paste, pattern.temp[,k])

  pattern.kinds <- names(table(pattern.paste))
  for(l in 1:length(pattern.kinds)){
    Cl <- which(pattern.paste == pattern.kinds[l])
    Rl <- sum(P.est[Cl] * P.alpha[Cl])
    Il <- sum(P.alpha[Cl])
    P.Xj.alpha[Cl] <- (Rl + 1e-10) / (Il + 2e-10)
  }
  P.Xj.alpha[which(P.Xj.alpha < 0.0001)] <- 0.0001
  P.Xj.alpha[which(P.Xj.alpha > 0.9999)] <- 0.9999

  return(P.Xj.alpha)
}
