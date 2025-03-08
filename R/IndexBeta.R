
#' Calculate \eqn{\beta}
#'
#' @description
#' The function is able to calculate the \eqn{\beta} index for all items after fitting \code{CDM} or directly.
#'
#' @param Y A required \eqn{N} × \eqn{I} matrix or \code{data.frame} consisting of the responses of \code{N} individuals
#'          to \eqn{N} × \eqn{I} items. Missing values need to be coded as \code{NA}.
#' @param Q A required binary \eqn{I} × \eqn{K} matrix containing the attributes not required or required 
#'          master the items. The \code{i}th row of the matrix is a binary indicator vector indicating which
#'          attributes are not required (coded by 0) and which attributes are required (coded by 1) to master
#'          item \eqn{i}.
#' @param CDM.obj An object of class \code{CDM.obj}. Can be \code{NULL}, but when it is not \code{NULL}, it enables
#'                rapid validation of the Q-matrix without the need for parameter estimation.
#'                @seealso \code{\link[Qval]{CDM}}.
#' @param model Type of model to be fitted; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"},
#'              \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#'
#' @details
#'  For item \eqn{i} with the q-vector of the \eqn{c}-th (\eqn{c = 1, 2, ..., 2^{K}}) type, the \eqn{\beta} index is computed as follows:
#' 
#' \deqn{
#'    \beta_{ic} = \sum_{l=1}^{2^K} \left| \frac{r_{li}}{n_l} P_{ic}(\boldsymbol{\alpha}_{l}) - 
#'                 \left(1 - \frac{r_{li}}{n_l}\right) \left[1 - P_{ic}(\boldsymbol{\alpha}_{l})\right] \right|
#'               = \sum_{l=1}^{2^K} \left| \frac{r_{li}}{n_l} - \left[1 - P_{ic}(\boldsymbol{\alpha}_{l}) \right] \right|
#'  }
#'  
#' In the formula, \eqn{r_{li}} represents the number of examinees in attribute mastery pattern \eqn{\boldsymbol{\alpha}_{l}} who correctly 
#' answered item \eqn{i}, while \eqn{n_l} is the total number of examinees in attribute mastery pattern \eqn{\boldsymbol{\alpha}_{l}}. 
#' \eqn{P_{ic}(\boldsymbol{\alpha}_{l})} denotes the probability that an examinee in attribute mastery pattern \eqn{\boldsymbol{\alpha}_{l}} answers 
#' item \eqn{i} correctly when the q-vector for item \eqn{i} is of the \eqn{c}-th type. In fact, 
#' \eqn{\frac{r_{li}}{n_l}} is the observed probability that an examinee in attribute mastery pattern \eqn{\boldsymbol{\alpha}_{l}} answers 
#' item \eqn{i} correctly, and \eqn{\beta_{jc}} represents the difference between the actual proportion of 
#' correct answers for item \eqn{i} in each attribute mastery pattern and the expected probability of answering the 
#' item incorrectly in that state. Therefore, to some extent, \eqn{\beta_{jc}} can be considered as a measure 
#' of discriminability.
#'
#'
#' @seealso \code{\link[Qval]{validation}}
#'
#' @return
#' An object of class \code{matrix}, which consisted of \eqn{\beta} index for each item and each possible attribute mastery pattern.
#'
#' @author
#' Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' Li, J., & Chen, P. (2024). A new Q-matrix validation method based on signal detection theory. British Journal of Mathematical and Statistical Psychology, 00, 1–33. DOI: 10.1111/bmsp.12371
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
#' 
#' model <- "DINA"
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ, model = model, distribute = "horder")
#'
#' ## calculate beta directly
#' beta <-get.beta(Y = example.data$dat, Q = example.Q, model = model)
#' print(beta)
#'
#' ## calculate beta after fitting CDM
#' example.CDM.obj <- CDM(example.data$dat, example.Q, model=model)
#' beta <-get.beta(CDM.obj = example.CDM.obj)
#' print(beta)
#'
#' @export
#' @importFrom GDINA attributepattern
#' 
#' 
#' @export
#' @importFrom GDINA attributepattern
#'
#'
get.beta <- function(Y = NULL, Q = NULL, CDM.obj = NULL, model = "GDINA"){
  
  if(is.null(CDM.obj) & (is.null(Y) | is.null(Q)))
    stop("one of [CDM.obj)] and [Y and Q] must not be NULL !!!")
  
  if(!is.null(Q)){
    I <- nrow(Q)
    K <- ncol(Q)
    L <- 2^K
    N <- nrow(Y)
  }else{
    I <- length(CDM.obj$analysis.obj$catprob.parm)
    K <- log2(length(CDM.obj$P.alpha))
    L <- length(CDM.obj$P.alpha)
    Y <- CDM.obj$analysis.obj$Y
    Q <- CDM.obj$analysis.obj$Q
    N <- nrow(Y)
  }
  
  # Create beta matrix
  beta <- matrix(NA, I, L)
  
  # Generate pattern names efficiently using apply and paste
  pattern <- attributepattern(K)
  pattern.names <- apply(pattern, 1, paste0, collapse = "")
  
  # Assign column and row names
  colnames(beta) <- pattern.names
  rownames(beta) <- paste0("item ", 1:I)
  
  if(is.null(CDM.obj))
    CDM.obj <- CDM(Y=Y, Q=Q, model=model)
  alpha.P <- CDM.obj$alpha.P
  P.alpha <- CDM.obj$P.alpha
  alpha <- CDM.obj$alpha
  P.alpha.Xi <- CDM.obj$P.alpha.Xi

  AMP <- as.matrix(personparm(CDM.obj$analysis.obj, what = "MAP")[, -(K+1)])
  beta_Ni_ri.obj <- beta_Ni_ri(pattern, AMP, Y)
  ri <- beta_Ni_ri.obj$ri
  Ni <- beta_Ni_ri.obj$Ni
  
  for(i in 1:I){
    P.est <- calculatePEst(Y[, i], P.alpha.Xi)
    P.Xi.alpha <- apply(pattern[-1, ], 1, function(row) P_GDINA(row, P.est, pattern, P.alpha))
    value <- abs((ri[, i] / Ni) * P.Xi.alpha - (1 - ri[, i] / Ni) * (1 - P.Xi.alpha))
    value[is.nan(value) | is.infinite(value)] <- 0
    beta[i, -1] <- colSums(value)
  }
  
  return(beta)
}