# set.seed(123)
# 
# library(Qval)
# library(ddpcr)
# 
# ## generate Q-matrix and data
# K <- 4
# I <- 20
# IQ <- list(
#   P0 = runif(I, 0.1, 0.3),
#   P1 = runif(I, 0.7, 0.9)
# )
# 
# QRR1 <- QRR2 <- 0
# times <- 100
# for(t in 1:times){
# 
#   cat("======================", t, "/", times, "======================\n")
# 
#   Q <- sim.Q(K, I)
#   data <- sim.data(Q = Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")
#   Y <- data$da
# 
#   ## simulate random mis-specifications
#   MQ <- sim.MQ(Q, 0.2)
# 
#   CDM.obj <- CDM(data$dat, MQ)
#   GDINA.obj <- CDM.obj$analysis.obj
# 
#   obj1 <- validation(data$dat, MQ, CDM.obj, method = "Wald", search.method = "PAA",  iter.level = "item", maxitr = 30, eps=0.95)
#   obj2 <- validation(data$dat, MQ, CDM.obj, method = "Wald", search.method = "stepwise", iter.level = "item", maxitr = 30, eps=0.95)
# 
#   QRR1 <- QRR1 + zQRR(Q, obj1$Q.sug)
#   QRR2 <- QRR2 + zQRR(Q, obj2$Q.sug)
#   cat("----------------------", t, "/", times, "----------------------\n")
#   cat(t, ": QRR.1 =", round(QRR1/t, 3), " QRR.2 =", round(QRR2/t, 3), "\n")
# }
