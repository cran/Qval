
#' @importFrom GDINA attributepattern
validation.MLR <- function(Y, Q, 
                         CDM.obj=NULL, method="EM", mono.constraint=TRUE, model="GDINA",
                         search.method="ESA", iter.level = "test", maxitr=20, 
                         verbose = TRUE){

  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  L <- 2^K
  pattern <- attributepattern(K)

  Q.MLR <- Q
  Q.pattern <- Q.pattern.ini <- apply(Q.MLR, 1, function(x) get_Pattern(x, pattern))

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
        # Initialize variables
        fit.index.i <- rep(Inf, L)
        q.possible <- integer(0)  # Use integer() for better performance
        fit.index.ini[i] <- Inf   # Default initialization for fit.index.ini
        
        # Loop to evaluate MLR for each pattern
        for(l in 2:L){
          MLR.il <- get.MLR(alpha.P, Y[, i], pattern[l, ])
          val.aic <- MLR.il$aic
          val.p <- MLR.il$p
          val.r <- MLR.il$r
          
          fit.index.i[l] <- val.aic  # Store AIC for the current pattern
          
          if(all(val.p <= 0.01) && all(val.r > 0)){
            q.possible <- c(q.possible, l)  # Append to q.possible if conditions are met
          }
          
          if(l == Q.pattern.ini[i]){
            fit.index.ini[i] <- val.aic  # Store initial AIC if matching
          }
        }
        
        # Select the best pattern based on the available options
        if(length(q.possible) > 1){
          best <- which.min(fit.index.i[q.possible])
          Q.pattern.cur[i] <- q.possible[best]
          fit.index.cur[i] <- fit.index.i[Q.pattern.cur[i]]
        } else if(length(q.possible) == 1){
          Q.pattern.cur[i] <- q.possible
          fit.index.cur[i] <- fit.index.i[q.possible]
        } else {
          best <- which.min(fit.index.i)
          Q.pattern.cur[i] <- best
          fit.index.cur[i] <- fit.index.i[best]
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
              q.possible.cur <- get_Pattern(Q.i.cur, pattern)

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
          possible.vector.cur <- get_Pattern(Q.i.cur, pattern)
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

    validating.items <- which(Q.pattern.ini != Q.pattern.cur)
    fit.index.delta <- abs(fit.index.ini - fit.index.cur)
    if(length(validating.items) > 0) {
      if(iter.level == "item"){
        if(sum(fit.index.delta) > 0.00010){
          validating.items <- which.max(fit.index.delta)
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
      cat(paste0('Iter  =', sprintf("%4d", iter), "/", sprintf("%4d", maxitr), ","),
          change, 'items have changed,',
          paste0("\u0394AIC=", formatC(sum(fit.index.delta[validating.items]), digits = 5, format = "f")), "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }

  return(list(Q.original = Q, Q.sug = Q.MLR, 
              process = Q.pattern, priority=priority, 
              iter = iter - 1))

}
