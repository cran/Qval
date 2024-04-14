correctQ.MLR <- function(Y, Q, CDM.obj=NULL, method="BM", mono.constraint=TRUE, model="GDINA",
                         search.method="ESA", maxitr=20, verbose = TRUE){

  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)

  Q.MLR <- Q
  Q.pattern.ini <- rep(2, I)
  for(i in 1:I)
    Q.pattern.ini[i] <- get.Pattern(Q.MLR[i, ], pattern)
  Q.pattern <- Q.pattern.ini

  iter <- 0
  while(iter < maxitr){
    iter <- iter + 1
    priority <- NULL

    fit.index.ini <- rep(Inf, I)
    if(iter != 1 | is.null(CDM.obj))
      CDM.obj <- CDM(Y, Q.MLR, method=method, mono.constraint=mono.constraint, model=model, verbose = 0)
    alpha.P <- CDM.obj$alpha.P
    alpha <- CDM.obj$alpha

    fit.index.cur <- rep(0, I)
    Q.pattern.cur <- rep(2, I)
    for(i in 1:I){

      MLR.ini <- get.MLR(alpha.P, Y[, i], Q.MLR[i, ])
      fit.index.ini[i] <- MLR.ini$aic
      if(all(MLR.ini$p <= 0.01) && all(MLR.ini$r > 0)){
        fit.index.i <- fit.index.ini[i]
        q.possible <- Q.pattern.ini[i]
      }else{
        fit.index.i <- Inf
        q.possible <- 2
      }

      ######################################## ESA ########################################
      if(search.method == "ESA"){
        fit.index.i <- rep(Inf, L)
        q.possible <- c()
        for(l in 2:L){
          MLR.il <- get.MLR(alpha.P, Y[, i], pattern[l, ])
          if(all(MLR.il$p <= 0.01) && all(MLR.il$r > 0)){
            q.possible <- c(q.possible, l)
          }
          fit.index.i[l] <- MLR.il$aic
          if(l == Q.pattern.ini[i])
            fit.index.ini[i] <- MLR.il$aic
        }
        if(length(q.possible) > 1){
          Q.pattern.cur[i] <- q.possible[which.min(fit.index.i[q.possible])]
          fit.index.cur[i] <- fit.index.i[Q.pattern.cur[i]]
        }else if(length(q.possible) == 1){
          Q.pattern.cur[i] <- q.possible
          fit.index.cur[i] <- fit.index.i[q.possible]
        }else{
          Q.pattern.cur[i] <- which.min(fit.index.i)
          fit.index.cur[i] <- min(fit.index.i)
        }
      }

      ######################################## SSA ########################################
      if(search.method == "SSA"){
        Q.i <- rep(0, K)
        for(k in 1:K){
          q.possible.k <- NULL
          fit.index.i.k <- Inf
          Q.i.k <- Q.i
          for(kk in 1:K){
            Q.i.cur <- Q.i
            if(Q.i[kk] == 0){
              Q.i.cur[kk] <- 1
              q.possible.cur <- get.Pattern(Q.i.cur, pattern)

              MLR.il <- get.MLR(alpha.P, Y[, i], Q.i.cur)
              if(all(MLR.il$p <= 0.01) && all(MLR.il$r > 0)){
                if(MLR.il$aic < fit.index.i.k){
                  Q.i.k <- Q.i.cur
                  q.possible.k <- q.possible.cur
                  fit.index.i.k <- MLR.il$aic
                }
              }
            }
          }
          if(!is.null(q.possible.k)){
            if(fit.index.i.k <= fit.index.i){
              Q.i <- Q.i.k
              q.possible <- q.possible.k
              fit.index.i <- fit.index.i.k
            }
          }
        }
        Q.pattern.cur[i] <- q.possible
        fit.index.cur[i] <- fit.index.i
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

        fit.index.i <- q.possible <- c()
        for(k in 1:1:search.length){
          Q.i.cur <- Q.i
          att.posi <- which.max(priority.temp)
          Q.i.cur[att.posi] <- 1
          possible.vector.cur <- get.Pattern(Q.i.cur, pattern)
          priority.temp[att.posi] <- -Inf

          MLR.il <- get.MLR(alpha.P, Y[, i], Q.i.cur)
          if(all(MLR.il$p <= 0.01) && all(MLR.il$r > 0)){
            q.possible <- c(q.possible, possible.vector.cur)
            fit.index.i <- c(fit.index.i, MLR.il$aic)
          }
          Q.i <- Q.i.cur
        }
        if(length(q.possible) > 1){
          Q.pattern.cur[i] <- q.possible[which.min(fit.index.i)]
          fit.index.cur[i] <- min(fit.index.i)
        }else if (length(q.possible) == 1){
          Q.pattern.cur[i] <- q.possible
          fit.index.cur[i] <- fit.index.i
        }else if(length(q.possible) == 0){
          Q.pattern.cur[i] <- q.possible <- which.max(priority.cur) + 1
          fit.index.cur[i] <- get.MLR(alpha.P, Y[, i], pattern[q.possible, ])$aic
        }
      }

    }

    Q.pattern <- rbind(Q.pattern, Q.pattern.cur)
    if(iter > 2)
      if(all(Q.pattern.cur == Q.pattern[nrow(Q.pattern) - 2, ]))
        break
    fit.index.dif <- fit.index.ini - fit.index.cur
    validating.items <- which(fit.index.dif != 0)

    change <- 0
    isbreak <- FALSE
    for(i in validating.items){
      Q.temp <- Q.MLR
      Q.temp[i, ] <- pattern[Q.pattern.cur[i], ]
      if(all(colSums(Q.temp) > 0)){
        Q.MLR[i, ] <- pattern[Q.pattern.cur[i], ]
        Q.pattern.ini[i] <- Q.pattern.cur[i]
        change <- change + 1
      }else{
        isbreak <- TRUE
      }
    }
    if(change < 1)
      break
    if(isbreak){
      Q.MLR <- Q.temp
      break
    }
    if(verbose){
      cat(paste0('Iter = ', iter, "/", maxitr, ","), change, 'items have changed', "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }

  return(list(Q.original = Q, Q.sug = Q.MLR, priority=priority, iter = iter - 1))

}
