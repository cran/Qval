#' Calculate McFadden pseudo-R2
#'
#' @description
#' The function is able to calculate the McFadden pseudo-\eqn{R^{2}} (\eqn{R^{2}}) for all items after
#' fitting \code{CDM} or directly.
#'
#' @param Y A required \eqn{N} × \eqn{I} matrix or \code{data.frame} consisting of the responses of \code{N} individuals
#'          to \eqn{N} × \eqn{I} items. Missing values need to be coded as \code{NA}.
#' @param Q A required binary \eqn{I} × \eqn{K} matrix containing the attributes not required or required 
#'          master the items. The \code{i}th row of the matrix is a binary indicator vector indicating which
#'          attributes are not required (coded by 0) and which attributes are required (coded by 1) to master
#'          item \eqn{i}.
#' @param CDM.obj An object of class \code{CDM.obj}. Can can be NULL, but when it is not NULL, it
#'                enables rapid verification of the Q-matrix without the need for parameter estimation.
#'                @seealso \code{\link[Qval]{CDM}}.
#' @param model Type of model to fit; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"},
#'              \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#'
#' @details
#' The McFadden pseudo-\eqn{R^{2}} (McFadden, 1974) serves as a definitive model-fit index,
#' quantifying the proportion of variance explained by the observed responses. Comparable to the
#' squared multiple-correlation coefficient in linear statistical models, this coefficient of
#' determination finds its application in logistic regression models. Specifically, in the context
#' of the CDM, where probabilities of accurate item responses are predicted for each examinee,
#' the McFadden pseudo-\eqn{R^{2}} provides a metric to assess the alignment between these predictions
#' and the actual responses observed. Its computation is straightforward, following the formula:
#' \deqn{
#'  R_{i}^{2} = 1 - \frac{\log(L_{im})}{\log(L_{i0})}
#' }
#' where \eqn{\log(L_{im})} is the log-likelihood of the model, and  \eqn{\log(L_{i0})} is the log-likelihood of
#' the null model. If there were \eqn{N} examinees taking a test comprising \eqn{I} items, then \eqn{\log(L_{im})}
#' would be computed as:
#' \deqn{
#'  \log(L_{im}) =
#'  \sum_{p}^{N} \log \sum_{l=1}^{2^{K^\ast}} \pi(\boldsymbol{\alpha}_{l}^{\ast} | \boldsymbol{X}_{p})
#'    P_{i}(\boldsymbol{\alpha}_{l}^{\ast})^{X_{pi}} \left[ 1-P_{i}(\boldsymbol{\alpha}_{l}^{\ast}) \right] ^{(1-X_{pi})}
#' }
#' where \eqn{\pi(\boldsymbol{\alpha}_{l}^{\ast} | \boldsymbol{X}_{p})} is the posterior probability of examinee \eqn{p} with attribute
#' mastery pattern \eqn{\boldsymbol{\alpha}_{l}^{\ast}} when their response vector is \eqn{\boldsymbol{X}_{p}}, and \eqn{X_{pi}} is
#' examinee \eqn{p}'s response to item \eqn{i}. Let \eqn{\bar{X}_{i}} be the average probability of correctly responding
#' to item \eqn{i} across all \eqn{N} examinees; then \eqn{\log(L_{i0})} could be computed as:
#' \deqn{
#'  \log(L_{i0}) =
#'  \sum_{p}^{N} \log {\bar{X}_{i}}^{X_{pi}} {(1-\bar{X}_{i})}^{(1-X_{pi})}
#' }
#'
#' @seealso \code{\link[Qval]{validation}}
#'
#' @return
#' An object of class \code{matrix}, which consisted of \eqn{R^{2}} for each item and each possible attribute mastery pattern.
#'
#' @author
#' Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' McFadden, D. (1974). Conditional logit analysis of qualitative choice behavior. In P. Zarembka (Ed.), Frontiers in economics (pp.105–142). Academic Press.
#'
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2021). Balancing ft and parsimony to improve Q-matrix validation. British Journal of Mathematical and Statistical Psychology, 74, 110–130. DOI: 10.1111/bmsp.12228.
#'
#' Qin, H., & Guo, L. (2023). Using machine learning to improve Q-matrix validation. Behavior Research Methods. DOI: 10.3758/s13428-023-02126-0.
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' ## generate Q-matrix and data
#' K <- 3
#' I <- 20
#' example.Q <- sim.Q(K, I)
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")
#'
#' ## calculate R2 directly
#' R2 <-get.R2(Y = example.data$dat, Q = example.Q)
#' print(R2)
#'
#' ## calculate R2 after fitting CDM
#' example.CDM.obj <- CDM(example.data$dat, example.Q, model="GDINA")
#' R2 <-get.R2(CDM.obj = example.CDM.obj)
#' print(R2)
#'
#' @export
#' @importFrom GDINA attributepattern
#'

get.R2 <- function(Y = NULL, Q = NULL, CDM.obj = NULL, model = "GDINA") {
  if (is.null(CDM.obj) & (is.null(Y) | is.null(Q)))
    stop("one of [CDM.obj)] and [Y and Q] must not be NULL !!!")
  
  # Initialize dimensions and variables
  if (!is.null(Q)) {
    I <- nrow(Q)
    K <- ncol(Q)
    L <- 2^K
    N <- nrow(Y)
  } else {
    I <- length(CDM.obj$analysis.obj$catprob.parm)
    K <- log2(length(CDM.obj$P.alpha))
    L <- length(CDM.obj$P.alpha)
    Y <- CDM.obj$analysis.obj$Y
    Q <- CDM.obj$analysis.obj$Q
    N <- nrow(Y)
  }
  
  pattern <- attributepattern(K)
  
  # Ensure CDM object is initialized
  if (is.null(CDM.obj)) {
    CDM.obj <- CDM(Y, Q, model)
  }
  P.alpha <- CDM.obj$P.alpha
  P.alpha.Xi <- CDM.obj$P.alpha.Xi
  
  # Initialize matrices
  L.X <- matrix(0, I, L)
  P_est <- matrix(0, I, L)  # Precompute all P.est values
  
  # Precompute log-likelihoods
  for (i in 1:I) {
    P.i <- mean(Y[, i])
    L.X[i, 1] <- sum(log((P.i^Y[, i]) * ((1 - P.i)^(1 - Y[, i]))))
    
    # Precompute the P.est values
    P_est[i, ] <- calculatePEst(Y[, i], P.alpha.Xi)
    
    # Efficiently compute L.X for each pattern (excluding the first column) without loop
    P.Xi.alpha <- sapply(2:L, function(l) P_GDINA(pattern[l, ], P_est[i, ], pattern, P.alpha)) # Apply P_GDINA to each pattern
    L.X[i, 2:L] <- apply(P.Xi.alpha, 2, function(P_Xi) log_likelihood_i(Y[, i], P_Xi, P.alpha.Xi))  # Apply log_likelihood_i to each column of P.Xi.alpha
  }
  
  # Calculate R2
  R2 <- 1 - L.X / L.X[, 1]
  
  # Assign column and row names
  pattern.names <- apply(pattern, 1, paste, collapse = "")
  colnames(R2) <- pattern.names
  rownames(R2) <- paste0("item ", 1:I)
  
  return(R2)
}

