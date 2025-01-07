
#' 
#' @importFrom GDINA attributepattern
#' 
get.cov <- function(Y.i, P.alpha.Xi, P.Xi.alpha.reduced, Q.i.cur, pattern){
  att.posi <- which(Q.i.cur > 0)
  pattern.reduced <- attributepattern(length(att.posi))
  pattern.2 <- matrix(pattern[, att.posi], nrow = nrow(pattern))

  pattern.2.order <- rep(1, nrow(pattern.2))
  for(l in 1:nrow(pattern.2))
    pattern.2.order[l] <- get_Pattern(pattern.2[l, ], pattern.reduced)

  parde <- NULL
  for(l in 1:nrow(pattern.reduced)){
    P.alpha.Xi.l <- rowSums(P.alpha.Xi[, which(pattern.2.order == l)])
    parde.l <- P.alpha.Xi.l*(Y.i - P.Xi.alpha.reduced[l])/(P.Xi.alpha.reduced[l]*(1-P.Xi.alpha.reduced[l]))
    parde <- cbind(parde, parde.l)
  }

  infor.matr <- t(parde) %*% parde
  cov <- tryCatch({
    solve(infor.matr)
  }, error = function(e) {
    return(NULL)
  })

  pattern.reduced.names <- c()
  for(l in 1:nrow(pattern.reduced)){
    names.cur <- pattern.reduced[l, 1]
    if(ncol(pattern.reduced) > 1)
      for(k in 2:ncol(pattern.reduced))
        names.cur <- paste0(names.cur, pattern.reduced[l, k])
    pattern.reduced.names <- c(pattern.reduced.names, names.cur)
  }
  if(!is.null(cov))
    colnames(cov) <- rownames(cov) <- pattern.reduced.names
  return(cov)
}
