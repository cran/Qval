Wald.test <- function(Q.i, Q.i.cur, P.Xi.alpha.reduced, cov.cur){
  att.dif <- which(Q.i != Q.i.cur)
  att.posi <- which(Q.i.cur > 0)

  pattern.reduced <- attributepattern(length(att.posi))
  P.Xi.alpha.reduced <- t(t(P.Xi.alpha.reduced))

  Rmatrix <- matrix(0, 2^sum(Q.i), 2^sum(Q.i.cur))
  compare.order <- rep(0, nrow(pattern.reduced))
  comparr.posi <- 1
  for(l in 1:(nrow(pattern.reduced) - 1)){
    if(compare.order[l] != 0)
      next
    for(ll in (l+1):nrow(pattern.reduced))
      if(all(pattern.reduced[l, -att.dif] == pattern.reduced[ll, -att.dif])){
        compare.order[l] <- ll
        compare.order[ll] <- l
        Rmatrix[comparr.posi, sort(c(l, ll))] <- c(1, -1)
        comparr.posi <- comparr.posi + 1
      }
  }

  Wald.statistic <- tryCatch({
    t(Rmatrix %*% P.Xi.alpha.reduced) %*% solve(Rmatrix %*% cov.cur %*% t(Rmatrix)) %*% (Rmatrix %*% P.Xi.alpha.reduced)
  }, error = function(e) {
    return(NULL)
  })
  if(is.null(Wald.statistic))
    p.value <- 0
  else
    p.value <- 1 - stats::pchisq(abs(Wald.statistic), df=length(P.Xi.alpha.reduced))
  return(list(Wald.statistic=Wald.statistic, p.value=p.value))
}
