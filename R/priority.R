get.priority <- function(Y = NULL, Q = NULL, CDM.obj = NULL, model="GDINA"){

  if(is.null(CDM.obj) & (is.null(Y) | is.null(Q)))
    stop("one of [CDM.obj] and [Y, Q] must not be NULL !!!")

  if(!is.null(Q)){
    I <- nrow(Q)
  }else{
    Y <- CDM.obj$analysis.obj$Y
    Q <- CDM.obj$analysis.obj$Q
    I <- length(CDM.obj$analysis.obj$catprob.parm)
  }
  priority <- NULL

  if(is.null(CDM.obj))
    CDM.obj <- CDM(Y, Q, model)
  alpha.P <- CDM.obj$alpha.P

  for(i in 1:I){
    priority.cur <- get.MLRlasso(alpha.P, Y[, i])
    if(all(priority.cur <= 0))
      priority.cur[which.max(priority.cur)] <- 1
    priority <- rbind(priority, priority.cur)
  }

  return(priority)
}
