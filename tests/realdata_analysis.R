# library(Qval)
# library(CDM)
# library(openxlsx)
# 
# iter.max <- 1
# 
# # ########################### PISA 2000 ##################################
# data_names <- "PISA2000"
# Q_original <- as.matrix(data.pisa00R.ct$q.matrix)[, -c(7, 8)]
# Y <- as.matrix(data.pisa00R.ct$data[, 5:30])
# Y[, c(7, 12, 13)][which(Y[, c(7, 12, 13)] < 2)] <- 0
# Y[, c(7, 12, 13)][which(Y[, c(7, 12, 13)] > 1)] <- 1
# 
# J <- nrow(Q_original)
# K <- ncol(Q_original)
# 
# CDM_obj <- CDM(Y, Q_original, model="GDINA", method="BM", mono.constr = TRUE, verbose=1)
# 
# GDI_ESA <- validation(Y, Q_original, CDM.obj=CDM_obj, model="GDINA", method="GDI", search.method="ESA", maxitr=iter.max)
# GDI_PAA <- validation(Y, Q_original, CDM.obj=CDM_obj, model="GDINA", method="GDI", search.method="PAA", maxitr=iter.max)
# 
# Hull_ESA <- validation(Y, Q_original, CDM.obj=CDM_obj, model="GDINA", method="Hull", search.method="ESA", maxitr=iter.max)
# Hull_PAA <- validation(Y, Q_original, CDM.obj=CDM_obj, model="GDINA", method="Hull", search.method="PAA", maxitr=iter.max)
# 
# MLR_ESA <- validation(Y, Q_original, CDM.obj=CDM_obj, model="GDINA", method="MLR-B", search.method="ESA", maxitr=iter.max)
# MLR_PAA <- validation(Y, Q_original, CDM.obj=CDM_obj, model="GDINA", method="MLR-B", search.method="PAA", maxitr=iter.max)
# 
# model <- c("GDINA")
# fit_names <- c("npar", "-2LL", "AIC", "BIC", "CAIC", "SABIC",
#                "M2", "df", "M2.P", "RMSEA2", "SRMSR", "rate")
# fit_original <-
#   fit_GDI_ESA <- fit_GDI_PAA <-
#   fit_Hull_ESA <- fit_Hull_PAA <-
#   fit_MLR_ESA <- fit_MLR_PAA <- rep(0, length(fit_names))
# 
# cat("model:", model, "\n")
# 
# cat("original\n")
# fit_original <- c(fit(Y, Q_original, model), sum(abs(Q_original-Q_original))/(J*K))
# 
# cat("GDI ESA\n")
# fit_GDI_ESA <- c(fit(Y, GDI_ESA$Q.sug, model), sum(abs(GDI_ESA$Q.sug-Q_original))/(J*K))
# cat("GDI PAA\n")
# fit_GDI_PAA <- c(fit(Y, GDI_PAA$Q.sug, model), sum(abs(GDI_PAA$Q.sug-Q_original))/(J*K))
# 
# cat("Hull ESA\n")
# fit_Hull_ESA <- c(fit(Y, Hull_ESA$Q.sug, model), sum(abs(Hull_ESA$Q.sug-Q_original))/(J*K))
# cat("Hull PAA\n")
# fit_Hull_PAA <- c(fit(Y, Hull_PAA$Q.sug, model), sum(abs(Hull_PAA$Q.sug-Q_original))/(J*K))
# 
# cat("MLR ESA\n")
# fit_MLR_ESA <- c(fit(Y, MLR_ESA$Q.sug, model), sum(abs(MLR_ESA$Q.sug-Q_original))/(J*K))
# cat("MLR PAA\n")
# fit_MLR_PAA <- c(fit(Y, MLR_PAA$Q.sug, model), sum(abs(MLR_PAA$Q.sug-Q_original))/(J*K))
# 
# fit_obj <- rbind(unlist(fit_original), unlist(fit_GDI_ESA), unlist(fit_GDI_PAA), unlist(fit_Hull_ESA),
#                  unlist(fit_Hull_PAA), unlist(fit_MLR_ESA), unlist(fit_MLR_PAA))
# colnames(fit_obj) <- fit_names
# rownames(fit_obj) <- c("original", "GDI-ESA", "GDI-PAA", "Hull-ESA", "Hull-PAA", "MLR-ESA", "MLR-PAA")
# print(round(fit_obj, 3))
# 
# sum(abs(MLR_PAA$Q.sug - Q_original))
# 
# Q_matrix <- list(Q_original=Q_original,
#                  GDI_ESA=GDI_ESA$Q.sug, GDI_PAA=GDI_PAA$Q.sug,
#                  Hull_ESA=Hull_ESA$Q.sug, Hull_PAA=Hull_PAA$Q.sug,
#                  MLR_ESA=MLR_ESA$Q.sug, MLR_PAA=MLR_PAA$Q.sug)
# QRR <- matrix(0, length(Q_matrix), length(Q_matrix))
# rownames(QRR) <- colnames(QRR) <- c("original", paste0("GDI-", c("ESA", "PAA")),
#                                     paste0("Hull-", c("ESA", "PAA")), paste0("MLR-", c("ESA", "PAA")))
# for(q in 1:length(Q_matrix)){
#   for(qq in 1:length(Q_matrix))
#     QRR[q, qq] <- zQRR(Q_matrix[[q]], Q_matrix[[qq]])
# }
# 
# print(fit_obj)
# 
# # wb <- createWorkbook()
# # addWorksheet(wb, "fit index")
# # addWorksheet(wb, "original")
# # addWorksheet(wb, "GDI-ESA")
# # addWorksheet(wb, "GDI-PAA")
# # addWorksheet(wb, "Hull-ESA")
# # addWorksheet(wb, "Hull-PAA")
# # addWorksheet(wb, "MLR-ESA")
# # addWorksheet(wb, "MLR-PAA")
# # addWorksheet(wb, "QRR")
# # writeData(wb, "fit index", as.table(fit_obj))
# # writeData(wb, "original", as.table(Q_original))
# # writeData(wb, "GDI-ESA", as.table(GDI_ESA$Q.sug))
# # writeData(wb, "GDI-PAA", as.table(GDI_PAA$Q.sug))
# # writeData(wb, "Hull-ESA", as.table(Hull_ESA$Q.sug))
# # writeData(wb, "Hull-PAA", as.table(Hull_PAA$Q.sug))
# # writeData(wb, "MLR-ESA", as.table(MLR_ESA$Q.sug))
# # writeData(wb, "MLR-PAA", as.table(MLR_PAA$Q.sug))
# # writeData(wb, "QRR", as.table(QRR))
# # saveWorkbook(wb, paste0("results_realdata_", data_names, ".xlsx"),overwrite = TRUE)
# 
