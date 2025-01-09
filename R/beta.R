#'
#' @importFrom GDINA attributepattern
#' @import parallel
#'
correctQ.beta <- function(Y, Q, 
                          CDM.obj=NULL, method="EM", mono.constraint=TRUE, model="GDINA",
                          search.method="beta", maxitr=1, iter.level="test", 
                          criter="AIC", 
                          verbose = TRUE){
  
  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  pattern <- attributepattern(K)
  L <- nrow(pattern)

  Q.beta <- Q
  Q.pattern.ini <- rep(2, I)
  for(i in 1:I)
    Q.pattern.ini[i] <- get_Pattern(Q.beta[i, ], pattern)
  Q.pattern <- Q.pattern.ini

  best.pos <- NULL
  mod0 <- NULL
  
  iter <- 0
  while(iter < maxitr){
    iter <- iter + 1
    priority <- NULL
    
    if(iter != 1 | is.null(CDM.obj))
      CDM.obj <- CDM(Y, Q.beta, method=method, mono.constraint=mono.constraint, model=model, verbose = 0)
    alpha.P <- CDM.obj$alpha.P
    P.alpha <- CDM.obj$P.alpha
    alpha <- CDM.obj$alpha
    P.alpha.Xi <- CDM.obj$P.alpha.Xi
    
    Q.pattern.cur <- rep(2, I)
    fit.index.pre <- fit.index.cur <- rep(0, I)

    AMP <- as.matrix(personparm(CDM.obj$analysis.obj, what = "MAP")[, -(K+1)])
    beta_Ni_ri.obj <- beta_Ni_ri(pattern, AMP, Y)
    ri <- beta_Ni_ri.obj$ri
    Ni <- beta_Ni_ri.obj$Ni

    QvalEnv <- new.env()
    assign("Y", Y, envir = QvalEnv)
    assign("P.alpha.Xi", P.alpha.Xi, envir = QvalEnv)
    assign("P.alpha", P.alpha, envir = QvalEnv)
    assign("pattern", pattern, envir = QvalEnv)
    assign("ri", ri, envir = QvalEnv)
    assign("Ni", Ni, envir = QvalEnv)
    assign("Q.pattern.ini", Q.pattern.ini, envir = QvalEnv)
    assign("model", model, envir = QvalEnv)
    assign("criter", criter, envir = QvalEnv)
    assign("search.method", search.method, envir = QvalEnv)
    assign("P_GDINA", P.GDINA, envir = QvalEnv)
    assign("Q.beta", Q.beta, envir = QvalEnv)
    assign("L", L, envir = QvalEnv)
    assign("K", K, envir = QvalEnv)
    assign("alpha.P", alpha.P, envir = QvalEnv)
    assign("get.MLRlasso", get.MLRlasso, envir = QvalEnv)
    assign("priority", priority, envir = QvalEnv)
    assign("parallel_iter", parallel_iter, envir = QvalEnv)
    assign("parLapply", parLapply, envir = QvalEnv)

    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, c("Y", "P.alpha.Xi", "P.alpha", "pattern", "ri", "Ni", "Q.pattern.ini", "model", "criter",
                        "search.method", "P_GDINA", "Q.beta", "L", "K", "alpha.P", "get.MLRlasso", "priority"),
                  envir = QvalEnv) 
    results <- parLapply(cl, 1:I, 
                         fun = get("parallel_iter", envir = QvalEnv), 
                         Y = get("Y", envir = QvalEnv),
                         P.alpha.Xi = get("P.alpha.Xi", envir = QvalEnv),
                         P.alpha = get("P.alpha", envir = QvalEnv),
                         pattern = get("pattern", envir = QvalEnv),
                         ri = get("ri", envir = QvalEnv),
                         Ni = get("Ni", envir = QvalEnv),
                         Q.pattern.ini = get("Q.pattern.ini", envir = QvalEnv),
                         model = get("model", envir = QvalEnv),
                         criter = get("criter", envir = QvalEnv),
                         search.method = get("search.method", envir = QvalEnv),
                         P_GDINA = get("P_GDINA", envir = QvalEnv),
                         Q.beta = get("Q.beta", envir = QvalEnv),
                         L = get("L", envir = QvalEnv),
                         K = get("K", envir = QvalEnv),
                         alpha.P = get("alpha.P", envir = QvalEnv),
                         get.MLRlasso = get("get.MLRlasso", envir = QvalEnv),
                         priority = get("priority", envir = QvalEnv))
    stopCluster(cl)
    
    for(i in 1:I){
      fit.index.pre[i] <- results[[i]]$fit.index.pre
      fit.index.cur[i] <- results[[i]]$fit.index.cur
      Q.pattern.cur[i] <- results[[i]]$Q.pattern.cur
      priority <- rbind(priority, results[[i]]$priority)
      best.pos <- rbind(best.pos, results[[i]]$best.pos.i)
    }
    
    validating.items <- which(Q.pattern.ini != Q.pattern.cur)
    fit.index.delta <- abs(fit.index.cur - fit.index.pre)
    if(length(validating.items) > 0) {
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
          Q.pattern.cur[i] <- best.pos[validating.items[vi], att.vi]
        }
      }else if(iter.level == "item"){
        if(max(fit.index.delta) > 0.00010){
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
      Q.temp <- Q.beta
      Q.temp[i, ] <- pattern[Q.pattern.cur[i], ]
      if(all(colSums(Q.temp) > 0)){
        Q.beta[i, ] <- pattern[Q.pattern.cur[i], ]
        Q.pattern.ini[i] <- Q.pattern.cur[i]
        change <- change + 1
      }else{
        isbreak <- TRUE
      }
    }
    if(change < 1)
      break
    if(isbreak){
      Q.beta <- Q.temp
      break
    }
    
    if(verbose){
      cat(paste0('Iter  =', sprintf("%4d", iter), "/", sprintf("%4d", maxitr), ","),
          change, 'items have changed,',
          paste0(paste0("\u0394", criter, "="), formatC(sum(fit.index.delta[validating.items]), digits = 5, format = "f")), "\n")
    }
  }
  if(search.method == "PAA"){
    rownames(priority) <- rownames(Q)
    colnames(priority) <- colnames(Q)
  }
  
  return(list(Q.original = Q, Q.sug = Q.beta, 
              process = Q.pattern, priority=priority, 
              iter = iter - 1))
}

#' 
#' @title A tool for the \eqn{\beta} Method
#'
#' @description
#' This function performs a single iteration of the \eqn{\beta} method for A item's validation. It is designed 
#' to be used in parallel computing environments to speed up the validation process of the \eqn{\beta} method. 
#' The function is a utility function for \code{\link[Qval]{validation}}, and it should not be called independently by the user.
#'
#' @param i Item number that need to be validated.
#' @param Y Observed data matrix for validation.
#' @param P.alpha.Xi Individual posterior
#' @param P.alpha Attribute prior weights.
#' @param pattern The attribute mastery pattern matrix.
#' @param ri A vector that contains the numbers of examinees in each knowledge state who correctly answered item \eqn{i}.
#' @param Ni A vector that contains the total numbers of examinees in each knowledge state.
#' @param Q.pattern.ini Initial pattern number for the model.
#' @param model Model object used for fitting (e.g., GDINA).
#' @param criter Fit criterion ("AIC", "BIC", "CAIC", or "SABIC").
#' @param search.method Search method for model selection ("beta", "ESA", "SSA", or "PAA").
#' @param P_GDINA Function to calculate probabilities for GDINA model.
#' @param Q.beta Q-matrix for validation.
#' @param L Number of latent pattern.
#' @param K Number of attributes.
#' @param alpha.P Individuals' marginal mastery probabilities matrix (Tu et al., 2022)
#' @param get.MLRlasso Function for Lasso regression with multiple linear regression.
#' @param priority Vector of priorities for PAA method search.
#'
#' @return
#' An object of class \code{validation} is a \code{list} containing the following components:
#' \item{fit.index.pre}{The previous fit index value after applying the selected search method.}
#' \item{fit.index.cur}{The current fit index value after applying the selected search method.}
#' \item{Q.pattern.cur}{The pattern that corresponds to the optimal model configuration for the current iteration.}
#' \item{priority}{The priority vector used in the PAA method, if applicable.}
#'
#' @export
parallel_iter <- function(i, Y, P.alpha.Xi, P.alpha, pattern, ri, Ni, 
                          Q.pattern.ini, model, criter, search.method, 
                          P_GDINA, Q.beta, L, K, alpha.P, get.MLRlasso, 
                          priority) {
  result <- list()
  best.pos.i <- rep(NA, K)
  fit.index.K <- rep(Inf, K)
  
  q.possible <- 2
  P.est <- calculatePEst(Y[, i], P.alpha.Xi)
  mod0 <- GDINA(Y, pattern[Q.pattern.ini, ], model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))
  
  best.pos.i[sum(pattern[Q.pattern.ini[i], ])] <- get_Pattern(pattern[Q.pattern.ini[i], ], pattern)
  if (criter == "AIC") {
    result$fit.index.pre <- mod0$testfit$AIC
    result$fit.index.cur <- mod0$testfit$AIC
  } else if (criter == "BIC") {
    result$fit.index.pre <- mod0$testfit$BIC
    result$fit.index.cur <- mod0$testfit$BIC
  } else if (criter == "CAIC") {
    result$fit.index.pre <- mod0$testfit$CAIC
    result$fit.index.cur <- mod0$testfit$CAIC
  } else if (criter == "SABIC") {
    result$fit.index.pre <- mod0$testfit$SABIC
    result$fit.index.cur <- mod0$testfit$SABIC
  }
  
  ################################### beta search ###################################
  if (search.method == "beta") {
    beta <- rep(0, K)
    pattern.single <- which(rowSums(pattern) == 1)
    for (k in 1:length(pattern.single)) {
      P.Xi.alpha.cur <- P_GDINA(pattern[pattern.single[k], ], P.est, pattern, P.alpha)
      value <- abs((ri[, i] / Ni) * P.Xi.alpha.cur - (1 - ri[, i] / Ni) * (1 - P.Xi.alpha.cur))
      value[is.nan(value) | is.infinite(value)] <- 0
      beta[k] <- sum(value)
    }
    best.pos.i[1] <- which.max(beta) + 1
    
    att.max <- which.max(beta)
    att.min <- which.min(beta)
    pattern.search <- which(pattern[, att.max] == 1 & pattern[, att.min] == 0 | rowSums(pattern) == K)
    
    fit.index.temp <- rep(0, length(pattern.search))
    for (l in 1:length(pattern.search)) {
      Q.temp <- Q.beta
      Q.temp[i, ] <- pattern[pattern.search[l], ]
      mod0 <- GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))
      
      if (criter == "AIC") {
        fit.index.temp[l] <- mod0$testfit$AIC
      } else if (criter == "BIC") {
        fit.index.temp[l] <- mod0$testfit$BIC
      } else if (criter == "CAIC") {
        fit.index.temp[l] <- mod0$testfit$CAIC
      } else if (criter == "SABIC") {
        fit.index.temp[l] <- mod0$testfit$SABIC
      }
      
      if(fit.index.K[sum(pattern[pattern.search[l], ])] > fit.index.temp[l]){
        best.pos.i[sum(pattern[pattern.search[l], ])] <- l
        fit.index.K[sum(pattern[pattern.search[l], ])] <- fit.index.temp[l]
      }
    }
    
    q.possible <- pattern.search[which.min(fit.index.temp)]
    result$fit.index.cur <- fit.index.temp[which.min(fit.index.temp)]
  }

  ######################################## ESA ########################################
  if (search.method == "ESA") {
    fit.index.i <- rep(Inf, L)
    beta.i <- rep(-Inf, L)
    for (l in 2:L) {
      Q.temp <- Q.beta
      Q.temp[i, ] <- pattern[l, ]
      mod0 <- GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))
      
      if (criter == "AIC") {
        fit.index.i[l] <- mod0$testfit$AIC
      } else if (criter == "BIC") {
        fit.index.i[l] <- mod0$testfit$BIC
      } else if (criter == "CAIC") {
        fit.index.i[l] <- mod0$testfit$CAIC
      } else if (criter == "SABIC") {
        fit.index.i[l] <- mod0$testfit$SABIC
      }
      
      if(fit.index.K[sum(pattern[l, ])] > fit.index.i[l]){
        best.pos.i[sum(pattern[l, ])] <- l
        fit.index.K[sum(pattern[l, ])] <- fit.index.i[l]
      }
    }
    
    q.possible <- which.min(fit.index.i)
    result$fit.index.cur <- fit.index.i[q.possible]
  }

  ######################################## SSA ########################################
  if (search.method == "SSA") {
    Q.i <- rep(0, K)
    fit.index.i <- Inf
    for (k in 1:K) {
      q.possible.k <- NULL
      fit.index.i.k <- fit.index.i.k.temp <- Inf
      beta.i.k <- -Inf
      Q.i.k <- Q.i
      for (kk in 1:K) {
        Q.i.cur <- Q.i
        if (Q.i[kk] == 0) {
          Q.i.cur[kk] <- 1
          q.possible.cur <- get_Pattern(Q.i.cur, pattern)
          
          Q.temp <- Q.beta
          Q.temp[i, ] <- Q.i.cur
          mod0 <- GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))
          
          if (criter == "AIC") {
            fit.index.i.k.temp <- mod0$testfit$AIC
          } else if (criter == "BIC") {
            fit.index.i.k.temp <- mod0$testfit$BIC
          } else if (criter == "CAIC") {
            fit.index.i.k.temp <- mod0$testfit$CAIC
          } else if (criter == "SABIC") {
            fit.index.i.k.temp <- mod0$testfit$SABIC
          }
          
          if (fit.index.i.k.temp < fit.index.i.k) {
            Q.i.k <- Q.i.cur
            q.possible.k <- q.possible.cur
            fit.index.i.k <- fit.index.i.k.temp
          }
        }
      }
      if (!is.null(q.possible.k)) {
        if (fit.index.i.k <= fit.index.i) {
          Q.i <- Q.i.k
          q.possible <- q.possible.k
          fit.index.i <- fit.index.i.k
          best.pos.i[sum(Q.i)] <- q.possible.k
        }
      }
    }
    result$fit.index.cur <- fit.index.i
  }

  ######################################## PAA ########################################
  if (search.method == "PAA") {
    priority.cur <- get.MLRlasso(alpha.P, Y[, i])
    if (all(priority.cur <= 0)) {
      priority.cur[which.max(priority.cur)] <- 1
    }
    result$priority <- priority.cur
    
    priority.temp <- priority.cur
    Q.i <- rep(0, K)
    fit.index.i.k <- fit.index.i.k.temp <- Inf
    beta.i <- 0
    P.Xi.alpha.all <- matrix(0, nrow = length(P.est), ncol = K)
    
    search.length <- sum(priority.cur > 0)
    for (k in 1:search.length) {
      att.posi <- which.max(priority.temp)
      Q.i[att.posi] <- 1
      q.possible.cur <- get_Pattern(Q.i, pattern)
      priority.temp[att.posi] <- -Inf
      
      Q.temp <- Q.beta
      Q.temp[i, ] <- Q.i
      mod0 <- GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))
      
      if (criter == "AIC") {
        fit.index.i.k.temp <- mod0$testfit$AIC
      } else if (criter == "BIC") {
        fit.index.i.k.temp <- mod0$testfit$BIC
      } else if (criter == "CAIC") {
        fit.index.i.k.temp <- mod0$testfit$CAIC
      } else if (criter == "SABIC") {
        fit.index.i.k.temp <- mod0$testfit$SABIC
      }
      best.pos.i[sum(Q.i)] <- q.possible.cur
      
      if (fit.index.i.k > fit.index.i.k.temp) {
        fit.index.i.k <- fit.index.i.k.temp
        q.possible <- q.possible.cur
      }
    }
    result$fit.index.cur <- fit.index.i.k
  }
  result$Q.pattern.cur <- q.possible
  result$best.pos.i <- best.pos.i
  
  return(result)
}


