correctQ.Hull <- function(Y, Q, CDM.obj=NULL, method="BM", mono.constraint=TRUE, model="GDINA",
                          search.method="ESA", maxitr=1, iter.level="test",
                          criter="PVAF", verbose = TRUE){

  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)

  Q.Hull <- Q
  Q.pattern.ini <- rep(2, I)
  for(i in 1:I)
    Q.pattern.ini[i] <- get.Pattern(Q.Hull[i, ], pattern)
  Q.pattern <- Q.pattern.ini

  iter <- 0
  while(iter < maxitr){
    iter <- iter + 1
    priority <- NULL

    if(iter != 1 | is.null(CDM.obj))
      CDM.obj <- CDM(Y, Q.Hull, method=method, mono.constraint=mono.constraint, model=model, verbose = 0)
    alpha.P <- CDM.obj$alpha.P
    P.alpha <- CDM.obj$P.alpha
    alpha <- CDM.obj$alpha
    P.alpha.Xi <- CDM.obj$P.alpha.Xi

    Q.pattern.cur <- rep(2, I)
    criterion.pre <- criterion.cur <- rep(0, I)
    for(i in 1:I){
        P.est <- (colSums(Y[, i] * P.alpha.Xi) + 1e-10) / (colSums(P.alpha.Xi) + 2e-10)
      if(criter == "PVAF"){
        P.mean <- sum(P.est * P.alpha)
        zeta2.i.K <- sum((P.est - P.mean)^2 * P.alpha)

        P.Xi.alpha <- P.GDINA(pattern[Q.pattern.ini[i], ], P.est, pattern, P.alpha)
        criterion.pre[i] <- criterion.cur[i] <- sum((P.Xi.alpha - P.mean)^2 * P.alpha) / zeta2.i.K
      }else if(criter == "R2"){
        L.X <- rep(0, L)
        L.Xi <- rep(0, N)
        Y.temp <- matrix(Y[, i], N, L, byrow = FALSE)
        P.i <- mean(Y[, i])

        L.X[1] <- sum(log((P.i^Y[, i])*((1-P.i)^(1-Y[, i]))))

        P.Xi.alpha <- P.GDINA(pattern[Q.pattern.ini[i], ], P.est, pattern, P.alpha)
        P.Xi.alpha.temp <- matrix(P.Xi.alpha, N, L, byrow = TRUE)
        L.Xi <- rowSums(P.alpha.Xi*(P.Xi.alpha.temp^Y.temp)*((1-P.Xi.alpha.temp)^(1-Y.temp)))
        L.X[Q.pattern.ini[i]] <- sum(log(L.Xi))

        criterion.pre[i] <- criterion.cur[i] <- 1- L.X[Q.pattern.ini[i]] / L.X[1]
      }
      ######################################## ESA ########################################
      if(search.method == "ESA"){
        if(criter == "PVAF"){
          zeta2 <- rep(-Inf, K)
          zeta2.cur <- 0
          pattern.zeta2 <- rep(2, K)
          for(l in 2:L){
            P.Xi.alpha <- P.GDINA(pattern[l, ], P.est, pattern, P.alpha)
            zeta2.cur <- sum((P.Xi.alpha - P.mean)^2 * P.alpha)
            if(zeta2.cur > zeta2[sum(pattern[l, ])]){
              zeta2[sum(pattern[l, ])] <- zeta2.cur
              pattern.zeta2[sum(pattern[l, ])] <- l
            }
          }
          criterion <- zeta2 / zeta2[K]
        }else if(criter == "R2"){
          P.est.temp <- matrix(P.est, N, L, byrow = TRUE)
          for(l in 2:L){
            P.Xi.alpha <- P.GDINA(pattern[l, ], P.est, pattern, P.alpha)
            P.Xi.alpha.temp <- matrix(P.Xi.alpha, N, L, byrow = TRUE)
            L.Xi <- rowSums(P.alpha.Xi*(P.Xi.alpha.temp^Y.temp)*((1-P.Xi.alpha.temp)^(1-Y.temp)))
            L.X[l] <- sum(log(L.Xi))
          }
          criterion <- 1- L.X / L.X[1]
        }

        Hull <- rep(-Inf, K-1)
        fjk <- rep(0, K+1)
        npk <- rep(0, K+1)
        for(k in 1:K){
          fjk[k+1] <- criterion[k]
          npk[k+1] <- 2^k
        }
        for(k in 1:(K-1)){
          pre <- (fjk[k+1]-fjk[k]) / (npk[k+1]- npk[k])
          nex <- (fjk[k+2]-fjk[k+1]) / (npk[k+2]- npk[k+1])
          Hull[k] <- pre / nex
        }
        Hull <- keep.convex(Hull, fjk, npk)
        if(!any(is.nan(Hull))){
          if(Q.pattern.cur[i] == 2 & all(Hull < 0)){
            Q.pattern.cur[i] <- L
          }else{
            Q.pattern.cur[i] <- pattern.zeta2[which.max(Hull)]
          }
        }
        criterion.cur[i] <- criterion[sum(pattern[Q.pattern.cur[i], ])]

        # {
        #   Hull[is.nan(Hull)] <- 0
        #   number.of.parameters <- c(0, npk[c(which(Hull>0)+1, K+1)])
        #   fit.index <- c(0, criterion[c(which(Hull>0)+1, K)])
        #   plot(number.of.parameters, fit.index, type = "o", main = "Hull plot")
        # }

      }

      ######################################## PAA ########################################
      if(search.method == "PAA"){
        priority.cur <- get.MLRlasso(alpha.P, Y[, i])
        if(all(priority.cur <= 0))
          priority.cur[which.max(priority.cur)] <- 1
        priority <- rbind(priority, priority.cur)

        priority.temp <- priority.cur
        Q.i <- rep(0, K)
        pattern.zeta2 <- c()

        search.length <- ifelse(length(which(priority.cur > 0)) < K, length(which(priority.cur > 0)) + 1, K)
        fjk <- npk <- Hull <- c(0)
        for(k in 1:search.length){
          Q.i.cur <- Q.i
          att.posi <- which.max(priority.temp)
          Q.i.cur[att.posi] <- 1
          q.possible.cur <- get.Pattern(Q.i.cur, pattern)
          priority.temp[att.posi] <- -Inf

          pattern.zeta2[sum(Q.i.cur)] <- get.Pattern(Q.i.cur, pattern)

          if(criter == "PVAF"){
            P.Xi.alpha.cur <- P.GDINA(Q.i.cur, P.est, pattern, P.alpha)
            zeta2.i.k.cur <- sum((P.Xi.alpha.cur - P.mean)^2 * P.alpha)
            fjk <- c(fjk, zeta2.i.k.cur/zeta2.i.K)
          }else if(criter == "R2"){
            P.Xi.alpha <- P.GDINA(Q.i.cur, P.est, pattern, P.alpha)
            P.Xi.alpha.temp <- matrix(P.Xi.alpha, N, L, byrow = TRUE)
            L.Xi <- rowSums(P.alpha.Xi*(P.Xi.alpha.temp^Y.temp)*((1-P.Xi.alpha.temp)^(1-Y.temp)))
            L.X[pattern.zeta2[sum(Q.i.cur)]] <- sum(log(L.Xi))
            fjk <- c(fjk, 1- L.X[pattern.zeta2[sum(Q.i.cur)]] / L.X[1])
          }
          npk <- c(npk, 2^sum(Q.i.cur))
          Q.i <- Q.i.cur
        }
        for(k in 1:(search.length-1)){
          pre <- (fjk[k+1]-fjk[k]) / (npk[k+1]- npk[k])
          nex <- (fjk[k+2]-fjk[k+1]) / (npk[k+2]- npk[k+1])
          Hull[k] <- pre / nex
        }

        if(length(pattern.zeta2) == 1){
          Q.pattern.cur[i] <- pattern.zeta2
        }else if(length(pattern.zeta2) == 2){
          if(max(Hull) > 1){
            Q.pattern.cur[i] <- pattern.zeta2[which.max(Hull)]
          }else{
            Q.pattern.cur[i] <- pattern.zeta2[2]
          }
        }else{
          Hull <- keep.convex(Hull, fjk, npk)
          if(!any(is.nan(Hull))){
            if(Q.pattern.cur[i] == 2 & all(Hull < 0)){
              Q.pattern.cur[i] <- L
            }else{
              Q.pattern.cur[i] <- pattern.zeta2[which.max(Hull)]
            }
          }
        }
        criterion.cur[i] <- fjk[which.max(Hull) + 1]
      }
    }

    Q.pattern <- rbind(Q.pattern, Q.pattern.cur)
    if(iter > 2)
      if(all(Q.pattern.cur == Q.pattern[nrow(Q.pattern) - 2, ]))
        break
    validating.items <- which(Q.pattern.ini != Q.pattern.cur)
    criterion.delta <- abs(criterion.pre - criterion.cur)
    if(iter.level == "item"){
      if(sum(criterion.delta) > 0.00010)
        validating.items <- which.max(criterion.delta)
      else
        validating.items <- integer(0)
    }

    change <- 0
    isbreak <- FALSE
    for(i in validating.items){
      Q.temp <- Q.Hull
      Q.temp[i, ] <- pattern[Q.pattern.cur[i], ]
      if(all(colSums(Q.temp) > 0)){
        Q.Hull[i, ] <- pattern[Q.pattern.cur[i], ]
        Q.pattern.ini[i] <- Q.pattern.cur[i]
        change <- change + 1
      }else{
        isbreak <- TRUE
      }
    }
    if(change < 1)
      break
    if(isbreak){
      Q.Hull <- Q.temp
      break
    }
    if(verbose){
      cat(paste0('Iter=', iter, "/", maxitr, ","),
          change, 'items have changed,',
          paste0(paste0("delta ", criter, "="), formatC(sum(criterion.delta), digits = 5, format = "f")), "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }

  return(list(Q.original = Q, Q.sug = Q.Hull, priority=priority, iter = iter - 1))

}
