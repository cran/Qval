#'
#' Priority of Attribute
#' 
#' @description
#' This function will provide the priorities of attributes for all items.
#' 
#' @details
#' The calculation of priorities is straightforward: the priority of an attribute is the 
#' regression coefficient obtained from a LASSO multinomial logistic regression, with the attribute 
#' as the independent variable and the response data from the subjects as the dependent variable.  
#' The formula is as follows:
#' 
#' \deqn{
#'  \log[\frac{P(X_{\pi} = 1 | \mathbf{\Lambda}_{p})}{P(X_{\pi} = 0 | \mathbf{\Lambda}_{p})}] = 
#'  logit[P(X_{\pi} = 1 | \mathbf{\Lambda}_{p})] = 
#'  \beta_{i0} + \beta_{i1} \Lambda_{p1} + \ldots + \beta_{ik} \Lambda_{pk} + \ldots + \beta_{iK} \Lambda_{pK}
#' }
#' 
#' The LASSO loss function can be expressed as:
#' 
#' \deqn{l_{lasso}(\mathbf{X}_i | \mathbf{\Lambda}) = l(\mathbf{X}_i | \mathbf{\Lambda}) - \lambda |\mathbf{\beta}_i|}
#' 
#' The priority for attribute \eqn{i} is defined as: \eqn{\mathbf{priority}_i = [\beta_{i1}, \ldots, \beta_{ik}, \ldots, \beta_{iK}]}
#' 
#' @param Y A required \code{N} × \code{I} matrix or data.frame consisting of the responses of \code{N} individuals
#'          to \code{I} items. Missing values need to be coded as \code{NA}.
#' @param Q A required binary \code{I} × \code{K} containing the attributes not required or required, 0 or 1,
#'            to master the items. The \code{i}th row of the matrix is a binary indicator vector indicating which
#'            attributes are not required (coded by 0) and which attributes are required (coded by 1) to
#'            master item \code{i}.
#' @param CDM.obj An object of class \code{CDM.obj}. When it is not NULL, it enables rapid verification
#'                of the Q-matrix without the need for parameter estimation. @seealso \code{\link[Qval]{CDM}}.
#' @param model Type of model to fit; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"}
#'              , \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#'              @seealso \code{\link[Qval]{CDM}}.
#'              
#' @returns A matrix containing all attribute priorities.
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
get.priority <- function(Y = NULL, Q = NULL, CDM.obj = NULL, model="GDINA"){

  if(is.null(CDM.obj) & (is.null(Y) | is.null(Q)))
    stop("one of [CDM.obj] and [Y, Q] must not be NULL !!!")

  if(!is.null(Q)){
    I <- nrow(Q)
  }else{
    Y <- CDM.obj$analysis.obj$Y
    Q <- CDM.obj$analysis.obj$Q
    I <- length(CDM.obj$analysis.obj$catprob.parm)
  }
  priority <- NULL

  if(is.null(CDM.obj))
    CDM.obj <- CDM(Y, Q, model)
  alpha.P <- CDM.obj$alpha.P

  for(i in 1:I){
    priority.cur <- get.MLRlasso(alpha.P, Y[, i])
    if(all(priority.cur <= 0))
      priority.cur[which.max(priority.cur)] <- 1
    priority <- rbind(priority, priority.cur)
  }
  rownames(priority) <- rownames(Q)
  colnames(priority) <- colnames(Q)
  

  return(priority)
}
