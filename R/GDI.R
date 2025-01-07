#'
#' @importFrom GDINA attributepattern
#'
correctQ.GDI <- function(Y, Q, 
                         CDM.obj=NULL, method="EM", mono.constraint=TRUE, model="GDINA",
                         search.method="ESA", maxitr=1, iter.level="test", eps=0.95,
                         verbose = TRUE){
  
  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)
  eps.value <- eps

  Q.GDI <- Q
  Q.pattern.ini <- rep(2, I)
  for(i in 1:I)
    Q.pattern.ini[i] <- get_Pattern(Q.GDI[i, ], pattern)
  Q.pattern <- Q.pattern.ini

  iter <- 0
  while(iter < maxitr){
    iter <- iter + 1
    priority <- NULL

    if(iter != 1 | is.null(CDM.obj))
      CDM.obj <- CDM(Y, Q.GDI, method=method, mono.constraint=mono.constraint, model=model, verbose = 0)
    alpha.P <- CDM.obj$alpha.P
    P.alpha <- CDM.obj$P.alpha
    alpha <- CDM.obj$alpha
    P.alpha.Xi <- CDM.obj$P.alpha.Xi

    Q.pattern.cur <- rep(2, I)
    PVAF.pre <- PVAF.cur <- rep(0, I)
    for(i in 1:I){

      P.est <- calculatePEst(Y[, i], P.alpha.Xi)
      P.mean <- sum(P.est * P.alpha)
      
      if(eps == "logit"){
        IQ <- 1 - P.est[1] - P.est[L]
        eps.eq <- -0.405 + 2.867*IQ + 4.840*10^4*N - 3.316*10^3*I
        eps.value <- exp(eps.eq) /(exp(eps.eq) + 1) 
      }
      
      P.Xi.alpha.L <- P_GDINA(rep(1, K), P.est, pattern, P.alpha)
      zeta2.i.K <- sum((P.Xi.alpha.L - P.mean)^2 * P.alpha)

      P.Xi.alpha <- P_GDINA(pattern[Q.pattern.ini[i], ], P.est, pattern, P.alpha)
      PVAF.pre[i] <- PVAF.cur[i] <- sum((P.Xi.alpha - P.mean)^2 * P.alpha) / zeta2.i.K

      ######################################## ESA ########################################
      if(search.method == "ESA"){
        zeta2 <- rep(-Inf, L)
        for(l in 2:L){
          P.Xi.alpha <- P_GDINA(pattern[l, ], P.est, pattern, P.alpha)
          zeta2[l] <- sum((P.Xi.alpha - P.mean)^2 * P.alpha)
        }
        PVAF <- zeta2 / zeta2[L]

        for(l in 2:L){
          if(PVAF[l] >= eps.value){
            Q.pattern.cur[i] <- l
            PVAF.cur[i] <- PVAF[l]
            break
          }
        }
      }

      ######################################## SSA ########################################
      if(search.method == "SSA"){
        Q.i <- rep(0, K)

        for(k in 1:K){
          q.possible.k <- NULL
          PVAF.i.k <- -Inf
          Q.i.k <- Q.i
          for(kk in 1:K){
            Q.i.cur <- Q.i
            if(Q.i[kk] == 0){
              Q.i.cur[kk] <- 1
              q.possible.cur <- get_Pattern(Q.i.cur, pattern)

              P.Xi.alpha.cur <- P_GDINA(Q.i.cur, P.est, pattern, P.alpha)
              zeta2.i.k.cur <- sum((P.Xi.alpha.cur - P.mean)^2 * P.alpha)
              if(PVAF.i.k < zeta2.i.k.cur/zeta2.i.K){
                Q.i.k <- Q.i.cur
                q.possible.k <- q.possible.cur
                PVAF.i.k <- zeta2.i.k.cur/zeta2.i.K
              }
            }
          }
          if(!is.null(q.possible.k)){
            Q.i <- Q.i.k
            q.possible <- q.possible.k
            if(PVAF.i.k >= eps.value){
              PVAF.cur[i] <- PVAF.i.k
              break
            }
          }
        }
        Q.pattern.cur[i] <- q.possible
      }

      ######################################## PAA ########################################
      if (search.method == "PAA") {
        priority.cur <- get.MLRlasso(alpha.P, Y[, i])
        if (all(priority.cur <= 0)) {
          priority.cur[which.max(priority.cur)] <- 1
        }
        priority <- rbind(priority, priority.cur)

        priority.temp <- priority.cur
        Q.i <- rep(0, K)
        PVAF.i.k <- -Inf
        P.Xi.alpha.all <- matrix(0, nrow = length(P.est), ncol = K)

        search.length <- sum(priority.cur > 0)

        for (k in 1:search.length) {
          att.posi <- which.max(priority.temp)
          Q.i[att.posi] <- 1
          q.possible.cur <- get_Pattern(Q.i, pattern)
          priority.temp[att.posi] <- -Inf

          P.Xi.alpha.cur <- P_GDINA(Q.i, P.est, pattern, P.alpha)
          zeta2.i.k.cur <- sum((P.Xi.alpha.cur - P.mean)^2 * P.alpha)

          if (PVAF.i.k < zeta2.i.k.cur/zeta2.i.K) {
            q.possible <- q.possible.cur
            PVAF.cur[i] <- PVAF.i.k <- zeta2.i.k.cur/zeta2.i.K
            if (PVAF.i.k >= eps.value) {
              break
            }
          }
        }

        Q.pattern.cur[i] <- q.possible
      }
    }

    
    validating.items <- which(Q.pattern.ini != Q.pattern.cur)
    PVAF.delta <- abs(PVAF.cur - PVAF.pre)
    if(iter.level == "item"){
      if(sum(PVAF.delta) > 0.00010){
        validating.items <- which.max(PVAF.delta)
        Q.pattern.cur[-validating.items] <- Q.pattern.ini[-validating.items]
        Q.pattern <- rbind(Q.pattern, Q.pattern.cur)
      }else{
        validating.items <- integer(0)
      }
    }else{
      Q.pattern <- rbind(Q.pattern, Q.pattern.cur)
    }

    change <- 0
    isbreak <- FALSE
    for(i in validating.items){
      Q.temp <- Q.GDI
      Q.temp[i, ] <- pattern[Q.pattern.cur[i], ]
      if(all(colSums(Q.temp) > 0)){
        Q.GDI[i, ] <- pattern[Q.pattern.cur[i], ]
        Q.pattern.ini[i] <- Q.pattern.cur[i]
        change <- change + 1
      }else{
        isbreak <- TRUE
      }
    }
    if(change < 1)
      break
    if(isbreak){
      Q.GDI <- Q.temp
      break
    }

    if(verbose){
      cat(paste0('Iter  =', sprintf("%4d", iter), "/", sprintf("%4d", maxitr), ","),
          change, 'items have changed,',
          paste0("\u0394PVAF=", formatC(sum(PVAF.delta[validating.items]), digits = 5, format = "f")), "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }

  return(list(Q.original = Q, Q.sug = Q.GDI, 
              process = Q.pattern, priority=priority, 
              iter = iter - 1))
}
