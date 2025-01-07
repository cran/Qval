# set.seed(123)
# 
# ## for Q-matrix validation procedures
# library(Qval)
# 
# ## for the dataset PISA2000
# library(CDM)
# 
# # ########################### PISA 2000 ##################################
# Q <- as.matrix(data.pisa00R.ct$q.matrix)[, -c(7, 8)]
# Y <- as.matrix(data.pisa00R.ct$data[, 5:30])
# Y[, c(7, 12, 13)][which(Y[, c(7, 12, 13)] < 2)] <- 0
# Y[, c(7, 12, 13)][which(Y[, c(7, 12, 13)] > 1)] <- 1
# 
# print(Q)
# 
# ######## GDI ########
# Qval.obj <- validation(Y, Q, model="GDINA",
#                        method="GDI", search.method="PAA",
#                        iter.level = "test", maxitr=1, verbose=1)
# print(Qval.obj$Q.sug)
# 
# ######## Hull ########
# Qval.obj <- validation(Y, Q, model="GDINA",
#                        method="Hull", search.method="ESA",
#                        iter.level = "item", maxitr=150, verbose=1)
# 
# head(Qval.obj$Q.sug)
# 
# plot(Qval.obj$Hull.fit, i=2)
# 
# ######## Wald ########
# Qval.obj <- validation(Y, Q, model="GDINA",
#                        method="Wald", alpha.level=0.05, search.method="stepwise",
#                        iter.level = "test", maxitr=1, verbose=1)
# 
# head(Qval.obj$Q.sug)
# 
# CDM.obj <- CDM(Y, Q)
# 
# q1 <- c(1, 0, 0, 1, 0, 0)
# q2 <- c(0, 0, 0, 1, 0, 0)
# Wald.test(CDM.obj, q1, q2, i=1)
# 
# ######## fit ########
# Qval.obj <- validation(Y, Q, model="GDINA",
#                        method="MLR-B", search.method="PAA",
#                        iter.level = "test", maxitr=1, verbose=1)
# 
# fit(Y, Q, model="GDINA")
# 
# ######## beta ########
# Qval.obj <- validation(Y, Q, model="GDINA",
#                        method="beta", search.method="PAA",
#                        criter="AIC",
#                        iter.level = "test", maxitr=1, verbose=1)
# 
# head(Qval.obj$Q.sug)
