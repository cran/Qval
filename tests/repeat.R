# set.seed(123)
# 
# library(Qval)
# library(ddpcr)
# 
# ## generate Q-matrix and data
# K <- 5
# I <- 20
# IQ <- list(
#   P0 = runif(I, 0.1, 0.3),
#   P1 = runif(I, 0.7, 0.9)
# )
# 
# QRR1 <- QRR2 <- 0
# times <- 1
# for(t in 1:times){
#   quiet({
#     Q <- sim.Q(K, I)
#     data <- sim.data(Q = Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")
#     Y <- data$da
#     
#     ## simulate random mis-specifications
#     MQ <- sim.MQ(Q, 0.1)
#     
#     CDM.obj <- CDM(data$dat, MQ)
#     GDINA.obj <- CDM.obj$analysis.obj
#     
#     obj1 <- validation(data$dat, MQ, CDM.obj, method = "Hull", search.method = "ESA") 
#     obj2 <- validation(data$dat, MQ, CDM.obj, method = "Hull", search.method = "PAA") 
#   })
#   
#   QRR1 <- QRR1 + zQRR(Q, obj1$Q.sug)
#   QRR2 <- QRR2 + zQRR(Q, obj2$Q.sug)
#   cat(t, ": QRR.1 =", round(QRR1/t, 3), " QRR.2 =", round(QRR2/t, 3), "\n")
# }
# 