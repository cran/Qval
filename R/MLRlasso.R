#'
#'
#' @importFrom glmnet glmnet
#' @importFrom Matrix drop0
#' @importFrom stats coef predict
#' 
get.MLRlasso <- function(alpha.P, Y.i) {
  fit <- glmnet(alpha.P, Y.i, family = "binomial", alpha = 1, nlambda = 50)

  lin.preds <- predict(fit, newx = alpha.P, type = "link", s = fit$lambda)
  predicted.probs <- 1 / (1 + exp(-lin.preds))
  Yi.matrix <- matrix(Y.i, nrow=length(Y.i), ncol=length(fit$lambda), byrow = FALSE)
  loglik <- colSums(Yi.matrix * log(predicted.probs) + (1 - Yi.matrix) * log(1 - predicted.probs))
  temp <- coef(fit) != 0
  k <- colSums(matrix(as.numeric(temp), nrow = nrow(temp), byrow = FALSE))
  bic.values <- -2 * loglik + log(length(Y.i)) * k
  best.coefs <- drop0(coef(fit, s = fit$lambda[which.min(bic.values)]))[-1, ]

  return(as.vector(best.coefs))
}

