#' Generate Response
#'
#' @description
#' randomly generate response matrix according to certain conditions,
#' including attributes distribution, item quality, sample size, Q-matrix and cognitive diagnosis models (CDMs).
#'
#' @param Q The Q-matrix. A random 30 × 5 Q-matrix (\code{\link[Qval]{sim.Q}}) will be used if \code{Q = NULL}.
#' @param N Sample size. Default = 500.
#' @param IQ A list containing two \eqn{I}-length vectors: \code{P0} and \code{P1}.  
#'           \code{P0} represents the probability of examinees who have not mastered any attributes  
#'           (\eqn{[00...0]}) correctly answering the item, while \code{P1} represents the probability  
#'           of examinees who have mastered all attributes (\eqn{[11...1]}) correctly answering the item.  
#' @param model Type of model to be fitted; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"},
#'              \code{"ACDM"}, \code{"LLM"}, or \code{"rRUM"}.
#' @param distribute Attribute distributions; can be \code{"uniform"} for the uniform distribution,
#'                   \code{"mvnorm"} for the multivariate normal distribution (Chiu, Douglas, & Li,
#'                   2009) and \code{"horder"} for the higher-order distribution (Tu et al., 2022).
#' @param control A list of control parameters with elements:
#' \itemize{
#'     \item \code{sigma}  A positive-definite symmetric matrix specifying the variance-covariance
#'                        matrix when \code{distribute = "mvnorm"}. Default = 0.5 (Chiu, Douglas, & Li, 2009).
#'     \item \code{cutoffs}  A vector giving the cutoff for each attribute when \code{distribute = "mvnorm"}.
#'                          Default = \eqn{k/(1+K)} (Chiu, Douglas, & Li, 2009).
#'     \item \code{theta} A vector of length N representing the higher-order ability for each examinee.
#'                       By default, generate randomly from the standard normal distribution (Tu et al, 2022).
#'     \item \code{a} The slopes for the higher-order model when \code{distribute = "horder"}.
#'                   Default = 1.5 (Tu et al, 2022).
#'     \item \code{b} The intercepts when \code{distribute = "horder"}. By default, select equally spaced
#'                   values between -1.5 and 1.5 according to the number of attributes (Tu et al, 2022).
#'  }
#' @param verbose Logical indicating to print information or not. Default is \code{TRUE}
#'
#' @return Object of class \code{sim.data}.
#' An \code{sim.data} object initially gained by \code{\link[GDINA]{simGDINA}} function form \code{GDINA} package.
#' Elements that can be extracted using method extract include:
#' \describe{
#'  \item{dat}{An \code{N} × \code{I} simulated item response matrix.}
#'  \item{Q}{The Q-matrix.}
#'  \item{attribute}{An \code{N} × \code{K} matrix for inviduals' attribute patterns.}
#'  \item{catprob.parm}{A list of non-zero success probabilities for each attribute mastery pattern.}
#'  \item{delta.parm}{A list of delta parameters.}
#'  \item{higher.order.parm}{Higher-order parameters.}
#'  \item{mvnorm.parm}{Multivariate normal distribution parameters.}
#'  \item{LCprob.parm}{A matrix of success probabilities for each attribute mastery pattern.}
#' }
#'
#' @author Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' Chiu, C.-Y., Douglas, J. A., & Li, X. (2009). Cluster Analysis for Cognitive Diagnosis: Theory and Applications. Psychometrika, 74(4), 633-665. DOI: 10.1007/s11336-009-9125-0.
#'
#' Tu, D., Chiu, J., Ma, W., Wang, D., Cai, Y., & Ouyang, X. (2022). A multiple logistic regression-based (MLR-B) Q-matrix validation method for cognitive diagnosis models:A confirmatory approach. Behavior Research Methods. DOI: 10.3758/s13428-022-01880-x.
#'
#' @examples
#'
#'################################################################
#'#                           Example 1                          #
#'#          generate data follow the uniform distrbution        #
#'################################################################
#' library(Qval)
#'
#' set.seed(123)
#'
#' K <- 5
#' I <- 10
#' Q <- sim.Q(K, I)
#'
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#'
#' data.obj <- sim.data(Q = Q, N = 100, IQ=IQ, model = "GDINA", distribute = "uniform")
#'
#' print(data.obj$dat)
#'
#'################################################################
#'#                           Example 2                          #
#'#          generate data follow the mvnorm distrbution         #
#'################################################################
#' set.seed(123)
#' K <- 5
#' I <- 10
#' Q <- sim.Q(K, I)
#'
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#'
#' cutoffs <- sample(qnorm(c(1:K)/(K+1)), ncol(Q))
#' data.obj <- sim.data(Q = Q, N = 10, IQ=IQ, model = "GDINA", distribute = "mvnorm",
#'                  control = list(sigma = 0.5, cutoffs = cutoffs))
#'
#' print(data.obj$dat)
#'
#'#################################################################
#'#                            Example 3                          #
#'#           generate data follow the horder distrbution         #
#'#################################################################
#' set.seed(123)
#' K <- 5
#' I <- 10
#' Q <- sim.Q(K, I)
#'
#' IQ <- list(
#'   P0 = runif(I, 0.0, 0.2),
#'   P1 = runif(I, 0.8, 1.0)
#' )
#'
#' theta <- rnorm(10, 0, 1)
#' b <- seq(-1.5,1.5,length.out=K)
#' data.obj <- sim.data(Q = Q, N = 10, IQ=IQ, model = "GDINA", distribute = "horder",
#'                  control = list(theta = theta, a = 1.5, b = b))
#'
#' print(data.obj$dat)
#'
#' @export
#' @importFrom GDINA attributepattern
#' @importFrom GDINA simGDINA
#' @importFrom stats runif qnorm rnorm sd
#' 
#'

sim.data <- function(Q=NULL, N=NULL, IQ=list(P0=NULL, P1=NULL),
                     model="GDINA", distribute="uniform", control = NULL,
                     verbose = TRUE){
  
  simCall <- match.call()

  if(is.null(Q))
    Q <- sim.Q(5, 30)
  K <- ncol(Q)
  I <- nrow(Q)
  if(is.null(N))
    N <- 500
  if(is.null(IQ$P0))
    IQ$P0 <- runif(I, 0, 0.3)
  if(is.null(IQ$P1))
    IQ$P1 <- runif(I, 0.7, 1)
  if(verbose){
    cat("distribute = ",distribute,"\n",
        "model = ",model,"\n",
        "number of attributes: ", K, "\n",
        "number of items: ", I, "\n",
        "num of examinees: ", N, "\n",
        "average of P0 = ", round(mean(IQ$P0), 3), "\n",
        "average of P1 = ", round(mean(IQ$P1), 3), "\n")
  }

  gs <- cbind(IQ$P0, 1 - IQ$P1)
  if(distribute == "mvnorm") {
    if(is.null(control$sigma)){
      sigma <- 0.5
    }else{
      sigma <- control$sigma
    }
    if(is.null(control$cutoffs)){
      cutoffs <- sample(qnorm(c(1:K)/(K+1)), ncol(Q), replace = FALSE)
    }else{
      cutoffs <-control$cutoffs
    }
    if(verbose){
      cat("sigma =", round(sigma, 3), "\n", "cutoffs =", round(cutoffs, 3), "\n")
    }

    vcov <- matrix(sigma,K,K)
    diag(vcov) <- 1
    data <- simGDINA(N, Q, gs.parm = gs, model = model, att.dist = "mvnorm",
                            gs.args = list(type = "random", mono.constraint = TRUE),
                            mvnorm.parm=list(mean = rep(0,K), sigma = vcov, cutoffs = cutoffs))
  }

  if(distribute == "horder") {
    if(is.null(control$theta)){
      theta <- rnorm(N, 0, 1)
    }else{
      theta <- control$theta
    }
    if(is.null(control$a)){
      a <- runif(K,1.5, 1.5)
    }else{
      a <- control$a
    }
    if(is.null(control$b)){
      b <- sample(seq(-1.5,1.5,length.out=K), K, replace = FALSE)
    }else{
      b <-control$b
    }
    if(verbose){
      cat("theta_mean = ", round(mean(theta), 3), ", theta_sd =", round(sd(theta), 3), "\n",
          "a = ", round(a, 3), "\n", "b = ", round(b, 3), "\n")
    }

    data <- simGDINA(N, Q, gs.parm = gs, model = model, att.dist = "higher.order",
                     gs.args = list(type = "random", mono.constraint = TRUE),
                     higher.order.parm = list(theta = theta,lambda = data.frame(a=a, b=b)))
  }

  if(distribute == "uniform")
    data <- simGDINA(N, Q, gs.parm = gs, model = model,
                            gs.args = list(type = "random", mono.constraint = TRUE))
  
   dat = data$dat
   Q = data$Q
   attribute = data$attribute
   catprob.parm = data$catprob.parm
   delta.parm = data$delta.parm
   higher.order.parm = data$higher.order.parm
   mvnorm.parm = data$mvnorm.parm
   LCprob.parm = data$LCprob.parm
   
   out = list(dat=dat, Q=Q, 
              attribute=attribute, 
              catprob.parm=catprob.parm, 
              delta.parm=delta.parm, 
              higher.order.parm=higher.order.parm, 
              mvnorm.parm=mvnorm.parm, 
              LCprob.parm=LCprob.parm, 
              call = simCall, 
              arguments = list(
                Q=Q, 
                N=N, 
                IQ=IQ,
                model=model, 
                distribute=distribute, 
                control = control,
                verbose = verbose
              ))
   
   class(out) <- "sim.data"
  
  return(out)
}
