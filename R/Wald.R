
#' @importFrom GDINA attributepattern Qval
correctQ.Wald <- function(Y, Q, CDM.obj=NULL, mono.constraint = TRUE, 
                          search.method="stepwise", iter.level = "test", maxitr=1,
                          eps=0.95, alpha.level=0.05, verbose = TRUE){

  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)

  Q.Wald <- Q
  Q.pattern.ini <- rep(2, I)
  for(i in 1:I)
    Q.pattern.ini[i] <- get_Pattern(Q.Wald[i, ], pattern)
  Q.pattern <- Q.pattern.ini

  iter <- 0
  while(iter < maxitr){
    iter <- iter + 1
    priority <- NULL

    if(iter != 1 | is.null(CDM.obj))
      CDM.obj <- CDM(Y, Q.Wald, mono.constraint=mono.constraint, "GDINA", verbose = 0)
    alpha.P <- CDM.obj$alpha.P
    P.alpha <- CDM.obj$P.alpha
    alpha <- CDM.obj$alpha
    P.alpha.Xi <- CDM.obj$P.alpha.Xi

    Q.pattern.cur <- rep(2, I)
    PVAF.pre <- PVAF.cur <- rep(0, I)
    
    for(i in 1:I){
      
      P.est <- calculatePEst(Y[, i], P.alpha.Xi)
      P.mean <- sum(P.est * P.alpha)
      P.Xi.alpha.L <- P_GDINA(rep(1, K), P.est, pattern, P.alpha)
      zeta2.i.K <- sum((P.Xi.alpha.L - P.mean)^2 * P.alpha)
      
      P.Xi.alpha <- P_GDINA(pattern[Q.pattern.ini[i], ], P.est, pattern, P.alpha)
      PVAF.pre[i] <- PVAF.cur[i] <- sum((P.Xi.alpha - P.mean)^2 * P.alpha) / zeta2.i.K
      
      ############################ stepwise or forward (SSA) #############################
      if(search.method == "stepwise" || search.method == "SSA" || search.method == "forward"){
        
        PVAF.K <- rep(0, K)
        Q.i <- rep(0, K)
        Q.i.full <- rep(1, K)
        for(k in 1:K){
          Q.i.cur <- Q.i
          Q.i.cur[k] <- 1
          P.Xj.alpha.cur <- P_GDINA(Q.i.cur, P.est, pattern, P.alpha)
          zeta2.i.k.cur <- sum((P.Xj.alpha.cur - P.mean)^2 * P.alpha)
          PVAF.K[k] <- zeta2.i.k.cur / zeta2.i.K
        }
        Q.i[which.max(PVAF.K)] <- 1
        PVAF.i <- max(PVAF.K)
        PVAF.cur[i] <- PVAF.i
        q.possible <- get_Pattern(Q.i, pattern)
        
        if(PVAF.i < eps){
          loop <- TRUE
          att.num <- sum(Q.i)
          while(loop && att.num < K){
            
            att.dif <- which(Q.i.full != Q.i)
            att.posi <- which(Q.i != 0)
            
            add.new <- remove.old <- NULL
            for(k in att.dif){
              Q.i.cur <- Q.i
              Q.i.cur[k] <- 1
              
              P.Xj.alpha.cur <- P_GDINA(Q.i.cur, P.est, pattern, P.alpha)
              zeta2.i.k.cur <- sum((P.Xj.alpha.cur - P.mean)^2 * P.alpha)
              PVAF.i.cur <- zeta2.i.k.cur / zeta2.i.K
              
              Wald.obj <- Wald.test(CDM.obj, Q.i, Q.i.cur, i=i)
              add.new <- rbind(add.new, c(k, Wald.obj$p.value, PVAF.i.cur))
              
              remove.old.cur <- NULL
              for(kk in att.posi){
                Q.i.temp <- Q.i.cur
                Q.i.temp[kk] <- 0
                Wald.obj <- Wald.test(CDM.obj, Q.i.temp, Q.i.cur, i=i)
                remove.old.cur <- c(remove.old.cur, Wald.obj$p.value)
              }
              remove.old <- rbind(remove.old, remove.old.cur)
            }
            # colnames(remove.old) <- att.posi
            # colnames(add.new) <- c("att", "p", "PVAF")
            
            att.operate <- which(add.new[, 2] < alpha.level)
            if(length(att.operate) > 0){
              add.remove <- cbind(add.new[att.operate, , drop=FALSE], remove.old[att.operate, , drop=FALSE])
              add.remove <- add.remove[which.max(add.remove[, 2]), ]
              
              Q.i[add.remove[1]] <- 1
              if(search.method == "stepwise"){
                temp <- add.remove[4:length(add.remove)]
                if(any(temp > alpha.level)){
                  att.remove <- as.numeric(names(temp))
                  Q.i[att.remove[which.max(temp)]] <- 0
                }
              }
              
              P.Xj.alpha.cur <- P_GDINA(Q.i, P.est, pattern, P.alpha)
              zeta2.i.k.cur <- sum((P.Xj.alpha.cur - P.mean)^2 * P.alpha)
              PVAF.i <- zeta2.i.k.cur / zeta2.i.K
              if(PVAF.i > eps){
                q.possible <- get_Pattern(Q.i, pattern)
                PVAF.cur[i] <- PVAF.i
                loop <- FALSE
              }
            }else{
              q.possible <- get_Pattern(Q.i, pattern)
              PVAF.cur[i] <- PVAF.i
              loop <- FALSE
            }
          }
        }
        Q.pattern.cur[i] <- q.possible
      }
      

      ######################################## PAA ########################################
      if(search.method == "PAA"){
        priority.cur <- get.MLRlasso(alpha.P, Y[, i])
        if(all(priority.cur <= 0))
          priority.cur[which.max(priority.cur)] <- 1
        priority <- rbind(priority, priority.cur)

        priority.temp <- priority.cur
        
        Q.i <- rep(0, K)
        search.length <- length(which(priority.cur > 0))
      
        for(k in 1:search.length){
          Q.i.cur <- Q.i
          att.posi <- which.max(priority.temp)
          Q.i.cur[att.posi] <- 1
          q.possible.cur <- get_Pattern(Q.i.cur, pattern)
          priority.temp[att.posi] <- -Inf
          P.Xj.alpha.cur <- P_GDINA(Q.i.cur, P.est, pattern, P.alpha)
          zeta2.i.k.cur <- sum((P.Xj.alpha.cur - P.mean)^2 * P.alpha)
          PVAF.i.k.cur <- zeta2.i.k.cur/zeta2.i.K

          if(PVAF.i.k.cur >= eps | search.length == 1){
            PVAF.cur[i] <- PVAF.i.k.cur
            Q.i <- Q.i.cur
            q.possible <- q.possible.cur
            break
          }
          if(k == 1){
            Q.i <- Q.i.cur
            q.possible <- q.possible.cur
            PVAF.cur[i] <- PVAF.i.k.cur
            next
          }
          
          Wald.obj <- Wald.test(CDM.obj, Q.i, Q.i.cur, i=i)
          if(Wald.obj$p.value < alpha.level){
            Q.i <- Q.i.cur
            q.possible <- q.possible.cur
            PVAF.cur[i] <- PVAF.i.k.cur
          }
        }
        Q.pattern.cur[i] <- q.possible
      }
      
    }

    validating.items <- which(Q.pattern.ini != Q.pattern.cur)
    PVAF.delta <- abs(PVAF.cur - PVAF.pre)
    if(length(validating.items) > 0){
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
    }
    
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
      cat(paste0('Iter  =', sprintf("%4d", iter), "/", sprintf("%4d", maxitr), ","),
          change, 'items have changed,',
          paste0("\u0394PVAF=", formatC(sum(PVAF.delta[validating.items]), digits = 5, format = "f")), "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }

  return(list(Q.original = Q, Q.sug = Q.Wald, priority=priority, iter = iter - 1))

}
