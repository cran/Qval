#'
#' @importFrom GDINA attributepattern
#'
validation.Hull <- function(Y, Q, 
                          CDM.obj=NULL, method="EM", mono.constraint=TRUE, model="GDINA",
                          search.method="ESA", maxitr=1, iter.level="test",
                          criter="PVAF", 
                          verbose = TRUE){

  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)

  Q.Hull <- Q
  Q.pattern <- Q.pattern.ini <- apply(Q.Hull, 1, function(x) get_Pattern(x, pattern))
  
  Hull.fit <- list()

  iter <- 0
  while(iter < maxitr){
    best.pos <- matrix(NA, I, K)
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
    pattern.criterion <- rep(2, K)
    for(i in 1:I){
      P.est <- calculatePEst(Y[, i], P.alpha.Xi)
      
      if(criter == "PVAF"){
        P.mean <- sum(P.est * P.alpha)
        P.Xi.alpha.L <- P_GDINA(rep(1, K), P.est, pattern, P.alpha)
        zeta2.i.K <- sum((P.Xi.alpha.L - P.mean)^2 * P.alpha)

        P.Xi.alpha <- P_GDINA(pattern[Q.pattern.ini[i], ], P.est, pattern, P.alpha)
        criterion.pre[i] <- criterion.cur[i] <- sum((P.Xi.alpha - P.mean)^2 * P.alpha) / zeta2.i.K
      }else if(criter == "R2"){
        L.X <- rep(0, L)
        L.Xi <- rep(0, N)
        P.i <- mean(Y[, i])

        L.X[1] <- sum(log((P.i^Y[, i])*((1-P.i)^(1-Y[, i]))))

        P.Xi.alpha <- P_GDINA(pattern[Q.pattern.ini[i], ], P.est, pattern, P.alpha)
        L.X[Q.pattern.ini[i]] <- log_likelihood_i(Y[, i], P.Xi.alpha, P.alpha.Xi)

        criterion.pre[i] <- criterion.cur[i] <- 1- L.X[Q.pattern.ini[i]] / L.X[1]
      }
      
      if(search.method == "ESA" | search.method == "SSA"){
        pattern.sum <- apply(pattern, 1, sum)
        
        ######################################## ESA ########################################
        if(search.method == "ESA"){
          if(criter == "PVAF"){
            zeta2 <- rep(-Inf, K)
            for(l in 2:L) {
              sum.pattern.l <- pattern.sum[l]
              
              P.Xi.alpha <- P_GDINA(pattern[l, ], P.est, pattern, P.alpha)
              zeta2.cur <- sum((P.Xi.alpha - P.mean)^2 * P.alpha)
              
              if(zeta2.cur > zeta2[sum.pattern.l]) {
                zeta2[sum.pattern.l] <- zeta2.cur
                best.pos[i, sum.pattern.l] <- pattern.criterion[sum.pattern.l] <- l
              }
            }
            criterion <- zeta2 / zeta2[K]
            
          }else if(criter == "R2"){
            P.est.temp <- matrix(P.est, N, L, byrow = TRUE)
            R2 <- rep(-Inf, K)
            R2.cur <- 0
            for(l in 2:L){
              sum.pattern.l <- pattern.sum[l]
              P.Xi.alpha <- P_GDINA(pattern[l, ], P.est, pattern, P.alpha)
              L.X[l] <- log_likelihood_i(Y[, i], P.Xi.alpha, P.alpha.Xi)
              
              R2.cur <- 1- L.X[l] / L.X[1]
              if(R2.cur > R2[sum.pattern.l]){
                R2[sum.pattern.l] <- R2.cur
                best.pos[i, sum.pattern.l] <- pattern.criterion[sum.pattern.l] <- l
              }
            }
            criterion <- R2
          }
        }
        
        ######################################## SSA ########################################
        if(search.method == "SSA"){
          Q.i <- criterion <- rep(0, K)
          for(k in 1:K){
            Q.i.k <- Q.i
            for(kk in 1:K){
              Q.i.cur <- Q.i
              if(Q.i[kk] == 0){
                Q.i.cur[kk] <- 1
                q.possible.cur <- get_Pattern(Q.i.cur, pattern)
                P.Xi.alpha.cur <- P_GDINA(Q.i.cur, P.est, pattern, P.alpha)
                
                if(criter == "PVAF"){
                  PVAF.i.k.cur <- sum((P.Xi.alpha.cur - P.mean)^2 * P.alpha) / zeta2.i.K
                  if(criterion[sum(Q.i.cur)] < PVAF.i.k.cur){
                    Q.i.k <- Q.i.cur
                    best.pos[i, sum(Q.i.cur)] <- pattern.criterion[sum(Q.i.cur)] <- q.possible.cur
                    criterion[sum(Q.i.cur)] <- PVAF.i.k.cur
                  }
                }else if(criter == "R2"){
                  L.X.cur <- log_likelihood_i(Y[, i], P.Xi.alpha.cur, P.alpha.Xi)
                  
                  R2.i.k.cur <- 1- L.X.cur / L.X[1]
                  if(criterion[sum(Q.i.cur)] < R2.i.k.cur){
                    Q.i.k <- Q.i.cur
                    best.pos[i, sum(Q.i.cur)] <- pattern.criterion[sum(Q.i.cur)] <- q.possible.cur
                    criterion[sum(Q.i.cur)] <- R2.i.k.cur
                  }
                }
              }
            }
            Q.i <- Q.i.k
          }
        }
        
        # Initialize the fjk and npk vectors, with the first elements as 0 and adjusted for criteria and powers of 2
        fjk <- c(0, criterion)  # Add 0 at the beginning of the criterion to align the indices
        npk <- c(0, 2^(1:K))    # Add 0 at the beginning, and create a sequence of 2^(1:K)
        
        # Calculate the first-order differences for fjk and npk using the diff() function
        Hull <- diff(fjk) / diff(npk)  # Compute the ratio of differences between successive elements of fjk and npk
        
        # Normalize the Hull values by dividing each value by the next one (element-wise division)
        Hull <- Hull[1:(K-1)] / Hull[2:K]  # Divide each element of Hull by the next one (index-wise)
        
        Hull.convex <- keep.convex(Hull, fjk, npk)
        if(!any(is.nan(Hull.convex$Hull))){
          if(Q.pattern.cur[i] == 2 & all(Hull.convex$Hull < 0)){
            Q.pattern.cur[i] <- L
          }else{
            Q.pattern.cur[i] <- pattern.criterion[which.max(Hull.convex$Hull)]
          }
        }
        criterion.cur[i] <- criterion[sum(pattern[Q.pattern.cur[i], ])]
        
        number.of.parameters <- c(0, 2^(1:K))
        fit.index <- c(0, criterion)
      }
      
      
      ######################################## PAA ########################################
      if(search.method == "PAA"){
        
        priority.cur <- get.MLRlasso(alpha.P, Y[, i])
        if(all(priority.cur <= 0))
          priority.cur[which.max(priority.cur)] <- 1
        priority <- rbind(priority, priority.cur)

        priority.temp <- priority.cur
        Q.i <- rep(0, K)
        pattern.criterion <- c()

        search.length <- ifelse(length(which(priority.cur > 0)) < K, length(which(priority.cur > 0)) + 1, K)
        fjk <- npk <- Hull <- c(0)
        for(k in 1:search.length){
          Q.i.cur <- Q.i
          att.posi <- which.max(priority.temp)
          Q.i.cur[att.posi] <- 1
          q.possible.cur <- get_Pattern(Q.i.cur, pattern)
          priority.temp[att.posi] <- -Inf

          best.pos[i, sum(Q.i.cur)] <- pattern.criterion[sum(Q.i.cur)] <- get_Pattern(Q.i.cur, pattern)

          if(criter == "PVAF"){
            P.Xi.alpha.cur <- P_GDINA(Q.i.cur, P.est, pattern, P.alpha)
            zeta2.i.k.cur <- sum((P.Xi.alpha.cur - P.mean)^2 * P.alpha)
            fjk <- c(fjk, zeta2.i.k.cur/zeta2.i.K)
          }else if(criter == "R2"){
            P.Xi.alpha <- P_GDINA(Q.i.cur, P.est, pattern, P.alpha)
            L.X[pattern.criterion[sum(Q.i.cur)]] <- log_likelihood_i(Y[, i], P.Xi.alpha, P.alpha.Xi)

            fjk <- c(fjk, 1- L.X[pattern.criterion[sum(Q.i.cur)]] / L.X[1])
          }
          npk <- c(npk, 2^sum(Q.i.cur))
          Q.i <- Q.i.cur
        }
        
        # Calculate the first-order differences of fjk and npk using diff()
        pre <- diff(fjk) / diff(npk)  # Calculate the difference ratio for fjk and npk (first-order)
        nex <- diff(fjk, lag = 2) / diff(npk, lag = 2)  # Calculate the second-order difference ratio for fjk and npk
        # Compute Hull by dividing the first-order difference by the second-order difference
        Hull <- pre[1:(search.length-1)] / nex  # Hull values are computed using pre and nex
        # Process Hull with keep.convex
        Hull.convex <- keep.convex(Hull, fjk, npk)  # Apply the keep.convex function to the Hull values
        
        if(!any(is.nan(Hull.convex$Hull))){
          if(all(Hull.convex$Hull < 0)){
            Q.pattern.cur[i] <- pattern.criterion[length(pattern.criterion)]
          }else{
            Q.pattern.cur[i] <- pattern.criterion[which.max(Hull.convex$Hull)]
          }
        }
        criterion.cur[i] <- fjk[which(pattern.criterion == Q.pattern.cur[i]) + 1]
        
        number.of.parameters <- c(0, 2^(1:length(pattern.criterion)))
        fit.index <- fjk
        
      }
      posi <- sort(unique(unlist(Hull.convex$posi)))
      Hull.fit[[i]] <- list(number.of.parameters=number.of.parameters, fit.index=fit.index, posi=posi, 
                            pattern.criterion=pattern.criterion, pattern=pattern, 
                            criter=criter, sug=criterion.cur[i])
    }

    validating.items <- which(Q.pattern.ini != Q.pattern.cur)
    criterion.delta <- abs(criterion.cur - criterion.pre)
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
        if(sum(criterion.delta) > 0.00010){
          validating.items <- which.max(criterion.delta)
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
      cat(paste0('Iter  =', sprintf("%4d", iter), "/", sprintf("%4d", maxitr), ","),
          change, 'items have changed,',
          paste0(paste0("\u0394", criter, "="), formatC(sum(criterion.delta[validating.items]), digits = 5, format = "f")), "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }

  return(list(Q.original = Q, Q.sug = Q.Hull, 
              process = Q.pattern, priority=priority, 
              Hull.fit = Hull.fit, 
              iter = iter - 1))

}
