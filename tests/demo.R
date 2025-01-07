# ################################################################
# #                           Example 1                          #
# #             The GDI method to validate Q-matrix              #
# ################################################################
# # set.seed(111)
# 
# library(Qval)
# 
# ## generate Q-matrix and data
# K <- 3
# I <- 20
# model <- "GDINA"
# Q <- example.Q <- sim.Q(K, I)
# IQ <- list(
#   P0 = runif(I, 0.0, 0.2),
#   P1 = runif(I, 0.8, 1.0)
# )
# example.data <- sim.data(Q = example.Q, N = 1000, IQ = IQ, model = model, distribute = "horder")
# Y <- example.data$dat
# 
# ## simulate random mis-specifications
# Q <- example.MQ <- sim.MQ(example.Q, 0.1)
# 
# ## using MMLE/EM to fit CDM model first
# example.CDM.obj <- CDM(example.data$dat, example.MQ, model = model)
# # GDINA.obj <- example.CDM.obj$analysis.obj
# #
# Q.beta.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "beta", search.method = "PAA", criter="AIC", iter.level = "item", maxitr = 10)
# # Q.beta.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "beta", search.method = "PAA", criter="BIC", iter.level = "item", maxitr = 20)
# print(zQRR(example.Q, Q.beta.obj$Q.sug))
# #
# # Q.GDI.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, model = model, method = "GDI", iter.level = "test", maxitr = 3)
# # Q.GDI.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, model = model, method = "GDI", iter.level = "item", maxitr = 150)
# # print(zQRR(example.Q, Q.GDI.obj$Q.sug))
# #
# # Q.Hull.obj <- validation(example.data$dat, Q, method = "Hull", search.method="ESA", iter.level = "test", maxitr = 1)
# # plot(Q.Hull.obj$Hull.fit, 19)
# # Q.Hull.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "Hull", search.method="PAA", iter.level = "item", maxitr = 150)
# # print(zQRR(example.Q, Q.Hull.obj$Q.sug))
# #
# # Q.Wald.obj <- validation(example.data$dat, example.MQ, method = "Wald", search.method="stepwise", iter.level = "item", maxitr = 150)
# # Q.Wald.obj <- validation(example.data$dat, example.MQ, method = "Wald", search.method="PAA", iter.level = "item", maxitr = 150)
# # print(zQRR(example.Q, Q.Wald.obj$Q.sug))
# #
# # Q.MLR.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "MLR-B", search.method="PAA", iter.level = "test", maxitr = 20)
# # Q.MLR.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "MLR-B", search.method="PAA", iter.level = "item", maxitr = 150)
# # print(zQRR(example.Q, Q.MLR.obj$Q.sug))
# 
