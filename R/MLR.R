get.MLR <- function(alpha.P, Y.i, q.vector) {

  att <- which(q.vector == 1)
  data <- data.frame(alpha.P[, att], Y = Y.i)
  names(data) <- c(paste0("A", att), "Y")

  formula <- paste("Y ~", paste(names(data)[-length(names(data))], collapse = " + "))

  MLR <- stats::glm(stats::as.formula(formula), data = data, family = "binomial")
  coeffs <- summary(MLR)$coefficients
  r <- coeffs[-1, 1]
  p <- coeffs[-1, 4]
  aic <- MLR$aic

  return(list(p = p, r = r, aic = aic))
}
