#'
#' @importFrom GDINA attributepattern
#'
validation.GDI <- function(Y, Q, 
                         CDM.obj=NULL, method="EM", mono.constraint=TRUE, model="GDINA",
                         search.method="ESA", maxitr=1, iter.level="test", eps=0.95,
                         verbose = TRUE){
  
  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)

  Q.GDI <- Q
  Q.pattern <- Q.pattern.ini <- apply(Q.GDI, 1, function(x) get_Pattern(x, pattern))

  iter <- 0
  while(iter < maxitr){
    best.pos <- matrix(NA, I, K)
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
        IQ <- P.est[1] - P.est[L]
        eps.eq <- -0.405 + 2.867*IQ + 4.840*10^(-4)*N - 3.316*10^(-3)*I
        eps.value <- exp(eps.eq) /(exp(eps.eq) + 1) 
      }else{
        eps.value <- eps
      }
      
      P.Xi.alpha.L <- P_GDINA(rep(1, K), P.est, pattern, P.alpha)
      zeta2.i.K <- sum((P.Xi.alpha.L - P.mean)^2 * P.alpha)

      P.Xi.alpha <- P_GDINA(pattern[Q.pattern.ini[i], ], P.est, pattern, P.alpha)
      PVAF.pre[i] <- PVAF.cur[i] <- sum((P.Xi.alpha - P.mean)^2 * P.alpha) / zeta2.i.K

      ######################################## ESA ########################################
      if(search.method == "ESA"){
        # Precompute zeta2 values for all patterns
        P.Xi.alpha <- sapply(2:L, function(l) P_GDINA(pattern[l, ], P.est, pattern, P.alpha))
        zeta2 <- c(-Inf, apply(P.Xi.alpha, 2, function(x) sum((x - P.mean)^2 * P.alpha)))
        PVAF <- zeta2 / zeta2[L]

        # Update PVAF.K and best.pos
        pattern.sum <- rowSums(pattern[1:L, ])  # Calculate sums of each pattern
        PVAF.K <- rep(-Inf, K)
        for(l in 2:L){
          pos <- pattern.sum[l]
          if(PVAF.K[pos] < PVAF[l]){
            PVAF.K[pos] <- PVAF[l]
            best.pos[i, pos] <- l
          }
        }

        # Find the first pattern where PVAF[l] >= eps.value
        q.possible <- best.pos[i, which(PVAF.K >= eps.value)]
        if (length(q.possible) > 0) {
          Q.pattern.cur[i] <- q.possible[1]
          PVAF.cur[i] <- PVAF[q.possible[1]]
        }
      }

      ######################################## SSA ########################################
      if(search.method == "SSA"){
        Q.i <- rep(0, K)
        already <- FALSE
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
            best.pos[i, sum(Q.i.k)] <- q.possible.k
            if(PVAF.i.k >= eps.value && !already){
              q.possible <- q.possible.k
              PVAF.cur[i] <- PVAF.i.k
              already <- TRUE
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
            best.pos[i, sum(Q.i)] <- q.possible.cur
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
    if(length(validating.items) > 0){
      if(iter.level == "test.att"){
        prov.Q <- pattern[Q.pattern.ini, ]
        cur.Q <- pattern[Q.pattern.cur, ]
        Ki.prov <- rowSums(prov.Q[validating.items, , drop = F])
        Ki.cur <- rowSums(cur.Q[validating.items, , drop = F])
        
        att.change <- ifelse(Ki.prov < Ki.cur, 1, ifelse(Ki.prov == Ki.cur, 0, -1))
        new.att <- Ki.prov + att.change
        for(vi in length(validating.items)){
          att.vi <- new.att[vi]
          while(is.na(best.pos[validating.items[vi], att.vi])) {
            att.vi <- att.vi - 1
          }
          Q.pattern.cur[validating.items[vi]] <- best.pos[validating.items[vi], att.vi]
        }
      }else if(iter.level == "item"){
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
      cat(paste0('Iter  =', sprintf("%3d", iter), "/", sprintf("%3d", maxitr), ","),
          sprintf("%3d", change), 'items have changed,',
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
