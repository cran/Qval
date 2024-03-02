#' Calculate data fit indeces
#'
#' @description
#' Calculate relative fit indices (-2LL, AIC, BIC, CAIC, SABIC) and absolute fit indices (\eqn{M_2} test)
#' using the \code{testfit} function in the \code{GDINA} package.
#'
#' @param Y A required \code{N} × \code{I} matrix or data.frame consisting of the responses of \code{N}
#'          individuals to \code{I} items. Missing values should be coded as NA.
#' @param Q A required binary \code{I} × \code{K} matrix containing the attributes not required or required
#'          , coded as 0 or 1, to master the items. The \code{i}th row of the matrix is a binary indicator vector
#'          indicating which attributes are not required (coded as 0) and which attributes are required
#'          (coded as 1) to master item \code{i}.
#' @param model Type of model to be fitted; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"},
#'              \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#'
#' @return
#' An object of class \code{list}. The list contains various fit indices:
#'
#' \item{npar}{The number of parameters.}
#' \item{-2LL}{The Deviance.}
#' \item{AIC}{The Akaike information criterion.}
#' \item{BIC}{The Bayesian information criterion.}
#' \item{CAIC}{The consistent Akaike information criterion.}
#' \item{SABIC}{The Sample size Adjusted BIC.}
#' \item{M2}{A vector consisting of \eqn{M_2} statistic, degrees of freedom, significance level, and \eqn{RMSEA_2} (Liu, Tian, & Xin, 2016).}
#' \item{SRMSR}{The standardized root mean squared residual (SRMSR; Ravand & Robitzsch, 2018).}
#'
#' @author
#' Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' Khaldi, R., Chiheb, R., & Afa, A.E. (2018). Feed-forward and Recurrent Neural Networks for Time Series Forecasting: Comparative Study. In: Proceedings of the International Conference on Learning and Optimization Algorithms: Theory and Applications (LOPAL 18). Association for Computing Machinery, New York, NY, USA, Article 18, 1–6. DOI: 10.1145/3230905.3230946.
#'
#' Liu, Y., Tian, W., & Xin, T. (2016). An application of M2 statistic to evaluate the fit of cognitive diagnostic models. Journal of Educational and Behavioral Statistics, 41, 3–26. DOI: 10.3102/1076998615621293.
#'
#' Ravand, H., & Robitzsch, A. (2018). Cognitive diagnostic model of best choice: a study of reading comprehension. Educational Psychology, 38, 1255–1277. DOI: 10.1080/01443410.2018.1489524.
#'
#' @examples
#' set.seed(123)
#'
#' library(Qval)
#'
#' ## generate Q-matrix and data to fit
#' K <- 5
#' I <- 30
#' example.Q <- sim.Q(K, I)
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")
#'
#' ## calculate fit indices
#' fit.indices <- fit(Y = example.data$dat, Q = example.Q, model = "GDINA")
#' print(fit.indices)
#'
#' @export
#' @importFrom GDINA modelfit
#' @importFrom GDINA npar
#'


fit <- function(Y, Q, model="GDINA"){
  GDINA.obj <- CDM(Y, Q, model = model, method = "EM")$analysis.obj
  testfit <- modelfit(GDINA.obj)
  if(is.null(testfit$M2)){
    testfit$M2 <- 0
    testfit$M2.pvalue <- 0
    testfit$M2.df <- 0
  }
  if(is.null(testfit$RMSEA2))
    testfit$RMSEA2 <- 0
  npar <- npar(GDINA.obj)[[1]]
  M2 <- c(testfit$M2, testfit$M2.df, testfit$M2.pvalue, testfit$RMSEA2)
  names(M2) <- c("M2", "df", "p.value", "RMSEA2")
  fit.index <- list(npar, testfit$Deviance, testfit$AIC, testfit$BIC, testfit$CAIC, testfit$SABIC,
                 M2, testfit$SRMSR)
  names(fit.index) <- c("npar", "-2LL", "AIC", "BIC", "CAIC", "SABIC",
                        "M2", "SRMSR")
  return(fit.index)
}
