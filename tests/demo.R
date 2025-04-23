# ################################################################
# #                           Example 1                          #
# #             The GDI method to validate Q-matrix              #
# ################################################################
# set.seed(321)
# 
# library(Qval)
# 
# ## generate Q-matrix and data
# K <- 5
# I <- 20
# model <- "GDINA"
# Q <- example.Q <- sim.Q(K, I)
# IQ <- list(
#   P0 = runif(I, 0.2, 0.4),
#   P1 = runif(I, 0.6, 0.8)
# )
# example.data <- sim.data(Q = example.Q, N = 1000, IQ = IQ, model = model, distribute = "horder")
# Y <- example.data$dat
# 
# # simulate random mis-specifications
# Q <- example.MQ <- sim.MQ(example.Q, 0.2)

# # using MMLE/EM to fit CDM model first
# CDM.obj <- example.CDM.obj <- CDM(example.data$dat, example.MQ, model = model)
# fit <- GDINA.obj <- example.CDM.obj$analysis.obj

# Q.beta.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "beta", search.method = "PAA", criter="AIC", iter.level = "test.att", maxitr = 1)
# Q.beta.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "beta", search.method = "beta", criter="BIC", iter.level = "test", maxitr = 20)
# print(zQRR(example.Q, Q.beta.obj$Q.sug))

# Q.GDI.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, model = model, method = "GDI", eps = "logit", search.method="ESA", iter.level = "test", maxitr = 1)
# Q.GDI.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, model = model, method = "GDI", search.method="ESA", iter.level = "test", maxitr = 1)
# Q.GDI.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, model = model, method = "GDI", iter.level = "item", maxitr = 150)
# print(zQRR(example.Q, Q.GDI.obj$Q.sug))

# Q.Hull.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "Hull", 
#                          criter = "PVAF", search.method="SSA", iter.level = "test", maxitr = 1)
# plot(Q.Hull.obj, 2)
# Q.Hull.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "Hull", search.method="PAA", iter.level = "item", maxitr = 150)
# print(zQRR(example.Q, Q.Hull.obj$Q.sug))

# Q.Wald.obj <- validation(example.data$dat, example.MQ, method = "Wald", search.method="SSA", iter.level = "test", maxitr = 1)
# Q.Wald.obj <- validation(example.data$dat, example.MQ, method = "Wald", search.method="PAA", iter.level = "item", maxitr = 150)
# print(zQRR(example.Q, Q.Wald.obj$Q.sug))

# Q.MLR.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "MLR-B", search.method="ESA", iter.level = "test", maxitr = 1)
# Q.MLR.obj <- validation(example.data$dat, example.MQ, example.CDM.obj, method = "MLR-B", search.method="PAA", iter.level = "item", maxitr = 150)
# print(zQRR(example.Q, Q.MLR.obj$Q.sug))

