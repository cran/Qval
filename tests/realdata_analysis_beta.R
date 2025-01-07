# library(Qval)
# 
# #---------PISA2000-----------
# 
# data(data.pisa00R.ct, package="CDM")
# dat <- data.pisa00R.ct$data 
# t_Q <- data.pisa00R.ct$q.matrix
# resp <- dat[, rownames(t_Q)]
# maxK <- apply( resp, 2, max, na.rm=TRUE )
# resp1 <- resp 
# for (ii in seq(1,ncol(resp)) ){ 
#   resp1[,ii] <- 1 * ( resp[,ii]==maxK[ii] ) 
# }
# 
# ORP <- resp1[,-c(7,9,12,13,14,26)]
# Q <- matrix(c(1,0,1,0,0, 1,0,1,1,0, 0,0,0,1,1, 0,1,1,1,0, 1,0,1,0,0, 
#               1,1,0,0,0, 1,1,1,0,0, 0,1,0,0,1, 0,1,1,0,0, 0,1,0,0,1,
#               1,1,0,0,1, 1,0,1,0,0, 1,0,0,0,0, 0,1,1,0,0, 0,1,0,0,0,
#               1,0,0,0,1, 0,1,1,0,0, 0,1,1,0,0, 0,1,0,1,0, 1,0,1,0,0),
#             20,5,byrow = TRUE)
# 
# #---------model selection-----------
# mod1 <- GDINA(dat = ORP, Q = Q, model = "DINA")
# mod2 <- GDINA(dat = ORP, Q = Q, model = "GDINA")
# anova(mod1, mod2)
# 
# #---------validation methods-----------
# model = "GDINA"
# GDINA.obj <- GDINA(ORP, Q, model, control = list(maxitr=300))
# Q.beta.obj <- validation(Y=ORP, Q=Q, method = "beta", search.method = "beta", criter="AIC", iter.level = "test", maxitr = 20)
# Q.Beta <- Q.beta.obj$Q.sug
# # Q.Beta              <- validateQ.Beta.gdina(GDINA.obj, ORP, Q, model)
# 
# #---------results of model fit-----------
# mod.Q        <- GDINA(ORP, Q, model, control = list(maxitr=300))
# mod.Beta     <- GDINA(ORP, Q.Beta, model, control = list(maxitr=300))
# 
# anova(mod.Q, mod.Beta)
# modelfit(mod.Q)
# modelfit(mod.Beta)
