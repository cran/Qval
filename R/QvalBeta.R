#'
#' @importFrom GDINA attributepattern
#' @import parallel
#'
validation.beta <- function(Y, Q, 
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
  Q.pattern <- Q.pattern.ini <- apply(Q.beta, 1, function(x) get_Pattern(x, pattern))
  
  mod0 <- NULL
  criter.index <- match(criter, c("AIC", "BIC", "CAIC", "SABIC"))
  
  iter <- 0
  while(iter < maxitr){
    iter <- iter + 1

    if(iter != 1 | is.null(CDM.obj)){
      CDM.obj <- CDM(Y, Q.beta, method=method, mono.constraint=mono.constraint, model=model, verbose = 0)
    }
    alpha.P <- CDM.obj$alpha.P
    P.alpha <- CDM.obj$P.alpha
    alpha <- CDM.obj$alpha
    P.alpha.Xi <- CDM.obj$P.alpha.Xi
    
    AMP <- as.matrix(personparm(CDM.obj$analysis.obj, what = "MAP")[, -(K+1)])
    beta_Ni_ri.obj <- beta_Ni_ri(pattern, AMP, Y)
    ri <- beta_Ni_ri.obj$ri
    Ni <- beta_Ni_ri.obj$Ni

    QvalEnv <- new.env()
    assign("Y", Y, envir = QvalEnv)
    assign("criter.index", criter.index, envir = QvalEnv)
    assign("P.alpha.Xi", P.alpha.Xi, envir = QvalEnv)
    assign("P.alpha", P.alpha, envir = QvalEnv)
    assign("pattern", pattern, envir = QvalEnv)
    assign("ri", ri, envir = QvalEnv)
    assign("Ni", Ni, envir = QvalEnv)
    assign("Q.pattern.ini", Q.pattern.ini, envir = QvalEnv)
    assign("model", model, envir = QvalEnv)
    assign("criter", criter, envir = QvalEnv)
    assign("search.method", search.method, envir = QvalEnv)
    assign("P_GDINA", P_GDINA, envir = QvalEnv)
    assign("Q.beta", Q.beta, envir = QvalEnv)
    assign("L", L, envir = QvalEnv)
    assign("K", K, envir = QvalEnv)
    assign("alpha.P", alpha.P, envir = QvalEnv)
    assign("get.MLRlasso", get.MLRlasso, envir = QvalEnv)
    assign("parallel_iter", parallel_iter, envir = QvalEnv)
    assign("parLapply", parLapply, envir = QvalEnv)

    cl <- makeCluster(detectCores() - 1)
    clusterExport(cl, c("Y", "P.alpha.Xi", "P.alpha", "pattern", "ri", "Ni", "Q.pattern.ini", "model", "criter",
                        "search.method", "P_GDINA", "Q.beta", "L", "K", "alpha.P", "get.MLRlasso"),
                  envir = QvalEnv) 
    results <- parLapply(cl, 1:I, 
                         fun = get("parallel_iter", envir = QvalEnv), 
                         Y = get("Y", envir = QvalEnv),
                         criter.index = get("criter.index", envir = QvalEnv),
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
                         get.MLRlasso = get("get.MLRlasso", envir = QvalEnv))
    stopCluster(cl)
    
    fit.index.pre <- sapply(results, function(x) x$fit.index.pre)
    fit.index.cur <- sapply(results, function(x) x$fit.index.cur)
    Q.pattern.cur <- sapply(results, function(x) x$Q.pattern.cur)
    priority <- do.call(rbind, lapply(results, function(x) x$priority))
    best.pos <- do.call(rbind, lapply(results, function(x) x$best.pos.i))
    
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
          Q.pattern.cur[validating.items[vi]] <- best.pos[validating.items[vi], att.vi]
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
#' This function performs a single iteration of the \eqn{\beta} method for one item's validation. It is designed 
#' to be used in parallel computing environments to speed up the validation process of the \eqn{\beta} method. 
#' The function is a utility function for \code{\link[Qval]{validation}}.
#' When the user calls the \code{\link[Qval]{validation}} function with \code{method = "beta"},  
#' \code{\link[Qval]{parallel_iter}} runs automatically, so there is no need for the user to call \code{\link[Qval]{parallel_iter}}.
#' It may seem that \code{\link[Qval]{parallel_iter}}, as an internal function, could better serve users.  
#' However, we found that the \code{Qval} package must export it to resolve variable environment conflicts in R  
#' and enable parallel computation. Perhaps a better solution will be found in the future.  
#'
#' @param i An integer indicating the item number that needs to be validated.
#' @param Y A matrix of observed data used for validation.
#' @param criter.index An integer representing the index of the criterion.
#' @param P.alpha.Xi A matrix representing individual posterior probability.
#' @param P.alpha A vector of attribute prior weights.
#' @param pattern A matrix representing the attribute mastery patterns.
#' @param ri A vector containing the number of examinees in each knowledge state who correctly answered item \eqn{i}.
#' @param Ni A vector containing the total number of examinees in each knowledge state.
#' @param Q.pattern.ini An integer representing the initial pattern order for the model.
#' @param model A model object used for fitting, such as the GDINA model.
#' @param criter A character string specifying the fit criterion. Possible values are "AIC", "BIC", "CAIC", or "SABIC".
#' @param search.method A character string specifying the search method for model selection. Options include "beta", "ESA", "SSA", or "PAA".
#' @param P_GDINA A function that calculates probabilities for the GDINA model.
#' @param Q.beta A Q-matrix used for validation.
#' @param L An integer representing the number of all attribute mastery patterns.
#' @param K An integer representing the number of attributes.
#' @param alpha.P A matrix of individuals' marginal mastery probabilities (Tu et al., 2022).
#' @param get.MLRlasso A function for Lasso regression with multiple linear regression.
#'
#' @return A \code{list} containing the following components:
#' \describe{
#'   \item{fit.index.pre}{The previous fit index value after applying the selected search method.}
#'   \item{fit.index.cur}{The current fit index value after applying the selected search method.}
#'   \item{Q.pattern.cur}{The pattern that corresponds to the optimal model configuration for the current iteration.}
#'   \item{priority}{The priority vector used in the PAA method, if applicable.}
#' }
#'
#' @importFrom GDINA GDINA
#' @export
parallel_iter <- function(i, Y, criter.index, P.alpha.Xi, P.alpha, pattern, ri, Ni, 
                          Q.pattern.ini, model, criter, search.method, 
                          P_GDINA, Q.beta, L, K, alpha.P, get.MLRlasso) {
  result <- list()
  best.pos.i <- rep(NA, K)
  fit.index.K <- rep(Inf, K)
  
  q.possible <- 2
  P.est <- calculatePEst(Y[, i], P.alpha.Xi)
  mod0 <- unlist(GDINA(Y, pattern[Q.pattern.ini, ], model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))$testfit)
  
  best.pos.i[sum(pattern[Q.pattern.ini[i], ])] <- get_Pattern(pattern[Q.pattern.ini[i], ], pattern)
  result$fit.index.pre <- result$fit.index.cur <- mod0[criter.index]
  
  ################################### beta search ###################################
  if (search.method == "beta") {
    ## determin the search space with beta
    pattern.single <- which(rowSums(pattern) == 1)
    P.Xi.alpha.cur <- apply(pattern[pattern.single, ], 1, function(row) P_GDINA(row, P.est, pattern, P.alpha))
    value <- abs((ri[, i] / Ni) * P.Xi.alpha.cur - (1 - ri[, i] / Ni) * (1 - P.Xi.alpha.cur))
    value[is.nan(value) | is.infinite(value)] <- 0
    beta <- colSums(value)
    best.pos.i[1] <- which.max(beta) + 1
    
    att.max <- which.max(beta)
    att.min <- which.min(beta)
    pattern.search <- which(pattern[, att.max] == 1 & pattern[, att.min] == 0 | rowSums(pattern) == K)

    ## search q-vectors in search space
    pattern.sum <- rowSums(pattern[pattern.search, ])
    fit.index.temp <- sapply(1:length(pattern.search), function(l) {
      Q.temp <- Q.beta
      Q.temp[i, ] <- pattern[pattern.search[l], ]
      mod0 <- unlist(GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))$testfit)
      
      if(fit.index.K[pattern.sum[l]] > mod0[criter.index]){
        best.pos.i[pattern.sum[l]] <<- pattern.search[l]
        fit.index.K[pattern.sum[l]] <<- mod0[criter.index]
      }
      
      return(mod0[criter.index])
    })

    q.possible <- pattern.search[which.min(fit.index.temp)]
    result$fit.index.cur <- fit.index.temp[which.min(fit.index.temp)]
  }

  ######################################## ESA ########################################
  if (search.method == "ESA") {
    pattern.sum <- rowSums(pattern[1:L, ])
    fit.index.i <- sapply(2:L, function(l) {
      Q.temp <- Q.beta
      Q.temp[i, ] <- pattern[l, ]
      mod0 <- unlist(GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))$testfit)
      
      if(fit.index.K[pattern.sum[l]] > mod0[criter.index]){
        best.pos.i[pattern.sum[l]] <<- l
        fit.index.K[pattern.sum[l]] <<- mod0[criter.index]
      }
      
      return(mod0[criter.index])
    })

    q.possible <- which.min(fit.index.i) + 1
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
          mod0 <- unlist(GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))$testfit)
          fit.index.i.k.temp <- mod0[criter.index]
          
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
      mod0 <- unlist(GDINA(Y, Q.temp, model, mono.constraint = TRUE, verbose=0, control = list(maxitr=300))$testfit)
      fit.index.i.k.temp <- mod0[criter.index]
      
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


