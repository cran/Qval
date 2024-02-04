
#' @importFrom GDINA attributepattern Qval
correctQ.Wald <- function(Y, Q, CDM.obj=NULL, model="GDINA",
                          search.method="ESA", maxitr=1,
                          eps=0.95, verbose = TRUE){

  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)

  Q.Wald <- Q
  Q.pattern.ini <- rep(2, I)
  for(i in 1:I)
    Q.pattern.ini[i] <- get.Pattern(Q.Wald[i, ], pattern)
  Q.pattern <- Q.pattern.ini

  iter <- 0
  while(iter < maxitr){
    iter <- iter + 1
    priority <- NULL

    if(iter != 1 | is.null(CDM.obj))
      CDM.obj <- CDM(Y, Q.Wald, "GDINA", verbose = 0)
    alpha.P <- CDM.obj$alpha.P
    P.alpha <- CDM.obj$P.alpha
    alpha <- CDM.obj$alpha
    P.alpha.Xi <- CDM.obj$P.alpha.Xi
    Mmatrix <- get.Mmatrix(pattern = pattern)

    Q.pattern.cur <- rep(2, I)
    ######################################## ESA ########################################
    if(search.method == "ESA"){
      Q.temp <- as.matrix(Qval(CDM.obj$analysis.obj, method = "wald")$sug.Q)
      for(i in 1:I){
        Q.pattern.cur[i] <- get.Pattern(Q.temp[i, ], pattern)
      }
    }

    ######################################## PAA ########################################
    if(search.method == "PAA"){
      for(i in 1:I){
        P.est <- (colSums(Y[, i] * P.alpha.Xi) + 1e-10) / (colSums(P.alpha.Xi) + 2e-10)
        P.mean <- sum(P.est * P.alpha)

        priority.cur <- get.MLRlasso(alpha.P, Y[, i])
        if(all(priority.cur <= 0))
          priority.cur[which.max(priority.cur)] <- 1
        priority <- rbind(priority, priority.cur)

        priority.temp <- priority.cur
        Q.i <- rep(0, K)
        zeta2.i.K <- sum((P.est - P.mean)^2 * P.alpha)
        PVAF.i.k <- -Inf
        q.possible <- which.max(priority.temp)+1

        search.length <- length(which(priority.cur > 0))

        for(k in 1:search.length){
          Q.i.cur <- Q.i
          att.posi <- which.max(priority.temp)
          Q.i.cur[att.posi] <- 1
          q.possible.cur <- get.Pattern(Q.i.cur, pattern)
          priority.temp[att.posi] <- -Inf

          P.Xj.alpha.cur <- P.GDINA(Q.i.cur, P.est, pattern, P.alpha)
          zeta2.i.k.cur <- sum((P.Xj.alpha.cur - P.mean)^2 * P.alpha)

          if(zeta2.i.k.cur/zeta2.i.K >= eps | search.length == 1){
            Q.i <- Q.i.cur
            q.possible <- q.possible.cur
            break
          }
          if(k == 1){
            Q.i <- Q.i.cur
            q.possible <- q.possible.cur
            next
          }
          P.Xi.alpha.reduced <- P.Xj.alpha.cur[which(Mmatrix[q.possible.cur, ] > 0)]
          cov.cur <- get.cov(Y[, i], P.alpha.Xi, P.Xi.alpha.reduced, Q.i.cur, pattern)
          if(!is.null(cov.cur)){
            Wald.obj <- Wald.test(Q.i, Q.i.cur, P.Xi.alpha.reduced, cov.cur)
            if(Wald.obj$p.value < 0.01){
              Q.i <- Q.i.cur
              q.possible <- q.possible.cur
            }
          }
        }
        Q.pattern.cur[i] <- q.possible
      }
    }

    Q.pattern <- rbind(Q.pattern, Q.pattern.cur)
    if(iter > 2)
      if(all(Q.pattern.cur == Q.pattern[nrow(Q.pattern) - 2, ]))
        break
    validating.items <- which(Q.pattern.ini != Q.pattern.cur)

    change <- 0
    isbreak <- FALSE
    for(i in validating.items){
      Q.temp <- Q.Wald
      Q.temp[i, ] <- pattern[Q.pattern.cur[i], ]
      if(all(colSums(Q.temp) > 0)){
        Q.Wald[i, ] <- pattern[Q.pattern.cur[i], ]
        Q.pattern.ini[i] <- Q.pattern.cur[i]
        change <- change + 1
      }else{
        isbreak <- TRUE
      }
    }
    if(change < 1)
      break
    if(isbreak){
      Q.Wald <- Q.temp
      break
    }
    if(verbose){
      cat(paste0('Iter = ', iter, "/", maxitr, ","), change, 'items changed', "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }

  return(list(Q.original = Q, Q.sug = Q.Wald, priority=priority, iter = iter - 1))

}
