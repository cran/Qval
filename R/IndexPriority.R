#'
#' Priority of Attribute
#' 
#' @description
#' This function will provide the priorities of attributes for all items.
#' 
#' @details
#' The calculation of priorities is straightforward (Qin & Guo, 2025): the priority of an attribute is the 
#' regression coefficient obtained from a LASSO multinomial logistic regression, with the attribute 
#' as the independent variable and the response data from the examinees as the dependent variable.  
#' The formula (Tu et al., 2022) is as follows:
#' 
#' \deqn{
#'  \log[\frac{P(X_{pi} = 1 | \boldsymbol{\Lambda}_{p})}{P(X_{pi} = 0 | \boldsymbol{\Lambda}_{p})}] = 
#'  logit[P(X_{pi} = 1 | \boldsymbol{\Lambda}_{p})] = 
#'  \beta_{i0} + \beta_{i1} \Lambda_{p1} + \ldots + \beta_{ik} \Lambda_{pk} + \ldots + \beta_{iK} \Lambda_{pK}
#' }
#' 
#' Where \eqn{X_{pi}} represents the response of examinee \eqn{p} on item \eqn{i},  
#' \eqn{\boldsymbol{\Lambda}_{p}} denotes the marginal mastery probabilities of examinee \eqn{p}  
#' (which can be obtained from the return value \code{alpha.P} of the \code{\link[Qval]{CDM}} function),  
#' \eqn{\beta_{i0}} is the intercept term, and \eqn{\beta_{ik}} represents the regression coefficient.  
#' 
#' The LASSO loss function can be expressed as:
#' 
#' \deqn{l_{lasso}(\boldsymbol{X}_i | \boldsymbol{\Lambda}) = l(\boldsymbol{X}_i | \boldsymbol{\Lambda}) - \lambda |\boldsymbol{\beta}_i|}
#' 
#' Where \eqn{l_{lasso}(\boldsymbol{X}_i | \boldsymbol{\Lambda})} is the penalized likelihood,  
#' \eqn{l(\boldsymbol{X}_i | \boldsymbol{\Lambda})} is the original likelihood,  
#' and \eqn{\lambda} is the tuning parameter for penalization (a larger value imposes a stronger penalty on 
#' \eqn{\boldsymbol{\beta}_i = [\beta_{i1}, \ldots, \beta_{ik}, \ldots, \beta_{iK}]}).  
#' The priority for attribute \eqn{i} is defined as: \eqn{\boldsymbol{priority}_i = \boldsymbol{\beta}_i = [\beta_{i1}, \ldots, \beta_{ik}, \ldots, \beta_{iK}]}
#' 
#' @param Y A required \eqn{N} × \eqn{I} matrix or \code{data.frame} consisting of the responses of \code{N} individuals
#'          to \eqn{N} × \eqn{I} items. Missing values need to be coded as \code{NA}.
#' @param Q A required binary \eqn{I} × \eqn{K} matrix containing the attributes not required or required 
#'          master the items. The \code{i}th row of the matrix is a binary indicator vector indicating which
#'          attributes are not required (coded by 0) and which attributes are required (coded by 1) to master
#'          item \eqn{i}.
#' @param CDM.obj An object of class \code{CDM.obj}. When it is not NULL, it enables rapid validation
#'                of the Q-matrix without the need for parameter estimation. @seealso \code{\link[Qval]{CDM}}.
#' @param model Type of model to fit; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"}
#'              , \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#'              @seealso \code{\link[Qval]{CDM}}.
#'              
#' @returns A matrix containing all attribute priorities.
#' 
#' @references
#' 
#' Qin, H., & Guo, L. (2025). Priority attribute algorithm for Q-matrix validation: A didactic. Behavior Research Methods, 57(1), 31. DOI: 10.3758/s13428-024-02547-5.
#'
#' Tu, D., Chiu, J., Ma, W., Wang, D., Cai, Y., & Ouyang, X. (2022). A multiple logistic regression-based (MLR-B) Q-matrix validation method for cognitive diagnosis models: A confirmatory approach. Behavior Research Methods. DOI: 10.3758/s13428-022-01880-x.
#'
#' @examples
#' set.seed(123)
#' library(Qval)
#' 
#' ## generate Q-matrix and data
#' K <- 5
#' I <- 20
#' IQ <- list(
#'   P0 = runif(I, 0.1, 0.3),
#'   P1 = runif(I, 0.7, 0.9)
#' )
#' 
#' \donttest{
#' Q <- sim.Q(K, I)
#' data <- sim.data(Q = Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")
#' MQ <- sim.MQ(Q, 0.1)
#' 
#' CDM.obj <- CDM(data$dat, MQ)
#' 
#' priority <- get.priority(data$dat, Q, CDM.obj)
#' head(priority)
#' }
#' 
#' 
#' @export
#'
get.priority <- function(Y = NULL, Q = NULL, CDM.obj = NULL, model="GDINA") {
  
  # Check input arguments
  if (is.null(CDM.obj) & (is.null(Y) | is.null(Q)))
    stop("One of [CDM.obj] and [Y, Q] must not be NULL !!!")
  
  # If Q is not NULL, extract related information; otherwise, extract from CDM.obj
  if (!is.null(Q)) {
    I <- nrow(Q)
  } else {
    Y <- CDM.obj$analysis.obj$Y
    Q <- CDM.obj$analysis.obj$Q
    I <- length(CDM.obj$analysis.obj$catprob.parm)
  }
  
  # If CDM.obj is NULL, initialize it
  if (is.null(CDM.obj)) {
    CDM.obj <- CDM(Y, Q, model)
  }
  alpha.P <- CDM.obj$alpha.P
  
  # Pre-allocate the priority matrix
  priority <- matrix(NA, nrow = I, ncol = ncol(Q))
  rownames(priority) <- rownames(Q)
  colnames(priority) <- colnames(Q)
  
  # Process each column of the matrix using a loop
  for (i in 1:I) {
    # Get the MLRlasso for the current column
    priority.cur <- get.MLRlasso(alpha.P, Y[, i])
    
    # If all values are non-positive, set the maximum to 1
    if (all(priority.cur <= 0)) {
      priority.cur[which.max(priority.cur)] <- 1
    }
    
    # Assign the result to the corresponding row in the priority matrix
    priority[i, ] <- priority.cur
  }
  
  return(priority)
}

