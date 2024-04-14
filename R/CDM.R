#' Parameter estimation for cognitive diagnosis models (CDMs) by MMLE/EM or MMLE/BM algorithm.
#'
#' @description
#' A function to estimate parameters for cognitive diagnosis models by MMLE/EM (de la Torre, 2009; de la Torre, 2011)
#' or MMLE/BM (Ma & Jiang, 2020) algorithm.The function imports various functions from the \code{GDINA} package,
#' parameter estimation for Cognitive Diagnostic Models was performed and extended. The \code{CDM} function not
#' only accomplishes parameter estimation for most commonly used models ( \code{GDINA}, \code{DINA}, \code{DINO},
#' \code{ACDM}, \code{LLM}, or \code{rRUM}) but also facilitates parameter estimation for the \code{LCDM}
#' model (Henson, Templin, & Willse, 2008; Tu et al., 2022). Furthermore, it incorporates Bayes modal estimation
#' (BM; Ma & Jiang, 2020) to obtain more reliable estimation results, especially in small sample sizes.
#' The monotonic constraints are able to be satisfied.
#'
#' @param Y A required \code{N} × \code{I} matrix or data.frame consisting of the responses of \code{N} individuals
#'          to × \code{I} items. Missing values need to be coded as \code{NA}.
#' @param Q A required binary \code{I} × \code{K} containing the attributes not required or required, 0 or 1, to
#'          master the items. The \code{i}th row of the matrix is a binary indicator vector indicating which
#'          attributes are not required (coded by 0) and which attributes are required (coded by 1) to master
#'          item \code{i}.
#' @param model Type of model to be fitted; can be \code{"GDINA"}, \code{"LCDM"}, \code{"DINA"}, \code{"DINO"}, \code{"ACDM"},
#'            \code{"LLM"}, or \code{"rRUM"}. Default = \code{"GDINA"}.
#' @param method  Type of mtehod to estimate CDMs' parameters; one out of \code{"EM"}, \code{"BM"}. Default = \code{"BM"}
#'                However, \code{"BM"} is only avaible when \code{method = "GDINA"}.
#' @param mono.constraint Logical indicating whether monotonicity constraints should be fulfilled in estimation.
#'                        Default = \code{TRUE}.
#' @param maxitr A vector for each item or nonzero category, or a scalar which will be used for all items
#'                 to specify the maximum number of EM or BM cycles allowed. Default = \code{2000}.
#' @param verbose Can be \code{0}, \code{1} or \code{2}, indicating to print no information, information
#'                for current iteration, or information for all iterations. Default = \code{1}.
#'
#' @details
#' CDMs are statistical models that fully integrates cognitive structure variables, which define the response
#' probability of subjects on questions by assuming the mechanism of action between attributes. In the
#' dichotomous test, this probability is the probability of answering correctly. According to the specificity
#' or generality of CDM assumptions, it can be divided into reduced CDM and saturated CDM.
#'
#' Reduced CDMs possess special and strong assumptions about the mechanisms of attribute interactions, leading
#' to clear interactions between attributes. Representative reduced models include the Deterministic Input,
#' Noisy and Gate (DINA) model (Haertel, 1989; Junker & Sijtsma, 2001; de la Torre & Douglas, 2004), the
#' Deterministic Input, Noisy or Gate (DINO) model (Templin & Henson, 2006), and the Additive Cognitive Diagnosis
#' Model (A-CDM; de la Torre, 2011), the reduced Reparametrized Unified Model (r-RUM; Hartz, 2002), among others.
#' Compared to reduced models, saturated models do not have strict assumptions about the mechanisms of attribute
#' interactions. When appropriate constraints are applied, they can be transformed into various reduced models
#' (Henson et al., 2008; de la Torre, 2011), such as the Log-Linear Cognitive Diagnosis Model (LCDM; Henson et
#' al., 2009) and the general Deterministic Input, Noisy and Gate model (G-DINA; de la Torre, 2011).
#'
#' The LCDM (Log-Linear Cognitive Diagnosis Model) is a saturated CDM fully proposed within the framework of
#' cognitive diagnosis. Unlike simplified models that only discuss the main effects of attributes, it also
#' considers the interactions between attributes, thus having more generalized assumptions about attributes.
#' Its definition of the probability of correct response is as follows:
#' \deqn{
#'    P(X_{pi}=1|\mathbf{\alpha}_{l}) =
#'    \frac{\exp(\lambda_{i0} + \mathbf{\lambda}_{i}^{T} \mathbf{h} (\mathbf{q_{i}}, \mathbf{\alpha_{l}}))}
#'    {1 + \exp(\lambda_{i0} + \mathbf{\lambda}_{i}^{T} \mathbf{h}(\mathbf{q_{i}}, \mathbf{\alpha_{l}}))}
#' }
#' \deqn{
#'    \mathbf{\lambda}_{i}^{T} \mathbf{h}(\mathbf{q_{i}}, \mathbf{\alpha_{l}}) =
#'    \lambda_{i0} + \sum_{k=1}^{K^\ast}\lambda_{ik}\alpha_{lk} +\sum_{k=1}^{K^\ast-1}\sum_{k'=k+1}^{K^\ast}
#'    \lambda_{ik}\lambda_{ik'}\alpha_{lk}\alpha_{lk'} +
#'    \cdots + \lambda_{12 \cdots K^\ast}\prod_{k=1}^{K^\ast}\alpha_{lk}
#' }
#' Where, \eqn{P(X_{pi}=1|\mathbf{\alpha}_{l})} represents the probability of a subject with attribute mastery
#' pattern \eqn{\mathbf{\alpha}_{l}}, where \eqn{l=1,2,\cdots,L} and \eqn{L=2^{K^\ast}}, correctly answering
#' item i.
#' Here, \eqn{K^\ast} denotes the number of attributes in the collapsed q-vector, \eqn{\lambda_{i0}} is the
#' intercept parameter, and \eqn{\mathbf{\lambda}_{i}=(\lambda_{i1}, \lambda_{i2}, \cdots, \lambda_{i12},
#' \cdots, \lambda_{i12{\cdots}K^\ast})} represents the effect vector of the attributes. Specifically,
#' \eqn{\lambda_{ik}} is the main effect of attribute \eqn{k}, \eqn{\lambda_{ikk'}} is the interaction effect between
#' attributes \eqn{k} and \eqn{k'}, and \eqn{\lambda_{j12{\cdots}K}} represents the interaction effect of all attributes.
#'
#' The general Deterministic Input, Noisy and Gate model (G-DINA), proposed by de la Torre (2011), is a saturated
#' model that offers three types of link functions: identity link, log link, and logit link, which are defined as follows:
#' \deqn{P(X_{pi}=1|\mathbf{\alpha}_{l}) =
#'    \delta_{i0} + \sum_{k=1}^{K^\ast}\delta_{ik}\alpha_{lk} +\sum_{k=1}^{K^\ast-1}\sum_{k'=k+1}^{K^\ast}\delta_{ik}\delta_{ik'}\alpha_{lk}\alpha_{lk'} +
#'    \cdots + \delta_{12{\cdots}K^\ast}\prod_{k=1}^{K^\ast}\alpha_{lk}
#' }
#' \deqn{log(P(X_{pi}=1|\mathbf{\alpha}_{l})) =
#'    v_{i0} + \sum_{k=1}^{K^\ast}v_{ik}\alpha_{lk} +\sum_{k=1}^{K^\ast-1}\sum_{k'=k+1}^{K^\ast}v_{ik}v_{ik'}\alpha_{lk}\alpha_{lk'} +
#'    \cdots + v_{12{\cdots}K^\ast}\prod_{k=1}^{K^\ast}\alpha_{lk}
#' }
#' \deqn{logit(P(X_{pi}=1|\mathbf{\alpha}_{l})) =
#'    \lambda_{i0} + \sum_{k=1}^{K^\ast}\lambda_{ik}\alpha_{lk} +\sum_{k=1}^{K^\ast-1}\sum_{k'=k+1}^{K^\ast}\lambda_{ik}\lambda_{ik'}\alpha_{lk}\alpha_{lk'} +
#'    \cdots + \lambda_{12{\cdots}K^\ast}\prod_{k=1}^{K^\ast}\alpha_{lk}
#' }
#' Where \eqn{\delta_{i0}}, \eqn{v_{i0}}, and \eqn{\lambda_{i0}} are the intercept parameters for the three
#' link functions, respectively; \eqn{\delta_{ik}}, \eqn{v_{ik}}, and \eqn{\lambda_{ik}} are the main effect
#' parameters of \eqn{\alpha_{lk}} for the three link functions, respectively; \eqn{\delta_{ikk'}}, \eqn{v_{ikk'}},
#' and \eqn{\lambda_{ikk'}} are the interaction effect parameters between \eqn{\alpha_{lk}} and \eqn{\alpha_{lk'}}
#' for the three link functions, respectively; and \eqn{\delta_{i12{\cdots }K^\ast}}, \eqn{v_{i12{\cdots}K^\ast}},
#' and \eqn{\lambda_{i12{\cdots}K^\ast}} are the interaction effect parameters of \eqn{\alpha_{l1}{\cdots}\alpha_{lK^\ast}}
#' for the three link functions, respectively. It can be observed that when the logit link is adopted, the
#' G-DINA model is equivalent to the LCDM model.
#'
#' Specifically, the A-CDM can be formulated as:
#' \deqn{P(X_{pi}=1|\mathbf{\alpha}_{l}) =
#'    \delta_{i0} + \sum_{k=1}^{K^\ast}\delta_{ik}\alpha_{lk}
#' }
#'
#' The RRUM, can be written as:
#' \deqn{log(P(X_{pi}=1|\mathbf{\alpha}_{l})) =
#'    \lambda_{i0} + \sum_{k=1}^{K^\ast}\lambda_{ik}\alpha_{lk}
#' }
#'
#' The item response function for LLM can be given by:
#' \deqn{logit(P(X_{pi}=1|\mathbf{\alpha}_{l})) =
#'    \lambda_{i0} + \sum_{k=1}^{K^\ast}\lambda_{ik}\alpha_{lk}
#' }
#'
#' In the DINA model, every item is characterized by two key parameters: guessing (g) and slip (s). Within
#' the traditional framework of DINA model parameterization, a latent variable \eqn{\eta}, specific to
#' individual \eqn{p} who has the attribute mastery pattern \eqn{\alpha_{l}} and item \eqn{i}, is defined as follows:
#' \deqn{
#'    \eta_{li}=\prod_{k=1}^{K}\alpha_{lk}^{q_{ik}}
#' }
#'
#' If individual \eqn{p} who has the attribute mastery pattern \eqn{\alpha_{l}} has acquired every attribute
#' required by item i, \eqn{\eta_{pi}} is given a value of 1. If not, \eqn{\eta_{pi}} is set to 0. The
#' DINA model's item response function can be concisely formulated as such:
#' \deqn{P(X_{pi}=1|\mathbf{\alpha}_{l}) =
#'    (1-s_j)^{\eta_{li}}g_j^{(1-\eta_{li})} =
#'    \delta_{i0}+\delta_{i12{\cdots}K}\prod_{k=1}^{K^\ast}\alpha_{lk}
#' }
#'
#' In contrast to the DINA model, the DINO model suggests that an individual can correctly respond to
#' an item if they have mastered at least one of the item's measured attributes. Additionally, like the
#' DINA model, the DINO model also accounts for parameters related to guessing and slipping. Therefore,
#' the main difference between DINO and DINA lies in their respective \eqn{\eta_{pi}} formulations. The
#' DINO model can be given by:
#' \deqn{\eta_{li} = 1-\prod_{k=1}^{K}(1 - \alpha_{lk})^{q_{lk}}}
#'
#' @return
#' An object of class \code{CDM.obj} is a \code{list} containing the following components:
#'
#' \item{analysis.obj}{An \code{GDINA} object gained from \code{GDINA} package or an \code{list} after BM algorithm,
#'                     depending on which estimation is used.}
#' \item{alpha}{Individuals' attribute parameters caculated by EAP method (Huebner & Wang, 2011)}
#' \item{P.alpha.Xi}{Individual posterior}
#' \item{alpha.P}{Individuals' marginal mastery probabilities matrix (Tu et al., 2022)}
#' \item{P.alpha}{Attribute prior weights for calculating marginalized likelihood in the last iteration}
#' \item{model.fit}{Some basic model-fit indeces, including \eqn{Deviance}, \eqn{npar}, \eqn{AIC}, \eqn{BIC}}
#'
#' @author Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' de la Torre, J. (2009). DINA Model and Parameter Estimation: A Didactic. Journal of Educational and Behavioral Statistics, 34(1), 115-130. DOI: 10.3102/1076998607309474.
#'
#' de la Torre, J., & Douglas, J. A. (2004). Higher-order latent trait models for cognitive diagnosis. Psychometrika, 69(3), 333-353. DOI: 10.1007/BF02295640.
#'
#' de la Torre, J. (2011). The Generalized DINA Model Framework. Psychometrika, 76(2), 179-199. DOI: 10.1007/s11336-011-9207-7.
#'
#' Haertel, E. H. (1989). Using restricted latent class models to map the skill structure of achievement items. Journal of Educational Measurement, 26(4), 301-323. DOI: 10.1111/j.1745-3984.1989.tb00336.x.
#'
#' Hartz, S. M. (2002). A Bayesian framework for the unified model for assessing cognitive abilities: Blending theory with practicality (Unpublished doctoral dissertation). University of Illinois at Urbana-Champaign.
#'
#' Henson, R. A., Templin, J. L., & Willse, J. T. (2008). Defining a Family of Cognitive Diagnosis Models Using Log-Linear Models with Latent Variables. Psychometrika, 74(2), 191-210. DOI: 10.1007/s11336-008-9089-5.
#'
#' Huebner, A., & Wang, C. (2011). A note on comparing examinee classification methods for cognitive diagnosis models. Educational and Psychological Measurement, 71, 407-419. DOI: 10.1177/0013164410388832.
#'
#' Junker, B. W., & Sijtsma, K. (2001). Cognitive assessment models with few assumptions, and connections with nonparametric item response theory. Applied Psychological Measurement, 25(3), 258-272. DOI: 10.1177/01466210122032064.
#'
#' Ma, W., & Jiang, Z. (2020). Estimating Cognitive Diagnosis Models in Small Samples: Bayes Modal Estimation and Monotonic Constraints. Applied Psychological Measurement, 45(2), 95-111. DOI: 10.1177/0146621620977681.
#'
#' Templin, J. L., & Henson, R. A. (2006). Measurement of psychological disorders using cognitive diagnosis models. Psychological methods, 11(3), 287-305. DOI: 10.1037/1082-989X.11.3.287.
#'
#' Tu, D., Chiu, J., Ma, W., Wang, D., Cai, Y., & Ouyang, X. (2022). A multiple logistic regression-based (MLR-B) Q-matrix validation method for cognitive diagnosis models: A confirmatory approach. Behavior Research Methods. DOI: 10.3758/s13428-022-01880-x.
#'
#' @seealso \code{\link[Qval]{validation}}.
#'
#' @examples
#'################################################################
#'#                           Example 1                          #
#'#            fit using MMLE/EM to fit the GDINA models         #
#'################################################################
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
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ,
#'                          model = "GDINA", distribute = "horder")
#'
#'\donttest{
#' ## using MMLE/EM to fit GDINA model
#' example.CDM.obj <- CDM(example.data$dat, example.Q, model = "GDINA",
#'                        method = "EM", maxitr = 2000, verbose = 1)
#'}
#'
#'
#'################################################################
#'#                           Example 2                          #
#'#               fit using MMLE/BM to fit the DINA              #
#'################################################################
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
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ,
#'                          model = "DINA", distribute = "horder")
#'
#'\donttest{
#' ## using MMLE/EM to fit GDINA model
#' example.CDM.obj <- CDM(example.data$dat, example.Q, model = "GDINA",
#'                        method = "BM", maxitr = 1000, verbose = 2)
#'}
#'
#'################################################################
#'#                           Example 3                          #
#'#              fit using MMLE/EM to fit the ACDM               #
#'################################################################
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
#' example.data <- sim.data(Q = example.Q, N = 500, IQ = IQ,
#'                          model = "ACDM", distribute = "horder")
#'
#'\donttest{
#' ## using MMLE/EM to fit GDINA model
#' example.CDM.obj <- CDM(example.data$dat, example.Q, model = "ACDM",
#'                        method = "EM", maxitr = 2000, verbose = 1)
#'}
#'
#' @export
#' @importFrom GDINA GDINA attributepattern designmatrix indlogPost LC2LG personparm
#'
#' @importFrom plyr llply
#' @importFrom nloptr slsqp
#' @importFrom Matrix drop0
#'

CDM <- function(Y, Q, model="GDINA", method="BM",
                mono.constraint=TRUE, maxitr=2000, verbose=1){

  if(all(method != c("EM", "BM")))
    stop("method must 'EM' or 'BM'!!! \n")

  if(method == "BM" & model != "GDINA")
    stop("model must be 'GDINA' when 'BM' !!! \n")

  Y <- as.matrix(Y)
  Q <- as.matrix(Q)
  N <- nrow(Y)
  I <- nrow(Q)
  K <- ncol(Q)
  KS <- attributepattern(K)
  L <- nrow(KS)

  if(method == "EM"){
    if(model=="LCDM"){
      mk <- rowSums(Q)
      numal <- 2^(mk)
      lcdm.ma <- list(designmatrix(Kj=mk[1], model = 'GDINA'))
      for (j in 2:nrow(Q)) {
        nextma <- designmatrix(Kj=mk[i], model = 'GDINA')
        list <- list(nextma)
        lcdm.ma <- c(lcdm.ma,list)
      }
      design.matrix2 <- lcdm.ma
      CDM.obj <- GDINA(Y, Q, verbose = verbose, mono.constraint = mono.constraint, control=list(maxitr = maxitr),
                              model = "UDF",linkfunc="logit", design.matrix = design.matrix2)
      if(verbose != 0)
        cat("\n")
    }else{
      CDM.obj <- GDINA(Y, Q, model = model, verbose = verbose, mono.constraint = mono.constraint, control=list(maxitr = maxitr))
      if(verbose != 0)
        cat("\n")
    }
    CDM.obj$Y <- Y
    CDM.obj$Q <- Q
    P.alpha.Xi <- exp(indlogPost(CDM.obj))
    P.alpha <- colSums(P.alpha.Xi) / sum(P.alpha.Xi)
    alpha.P <- personparm(CDM.obj, what = "mp")
    alpha <- KS[apply(P.alpha.Xi, 1, which.max), ]
  }
  if(method == "BM"){
    options('nloptr.show.inequality.warning'=FALSE)
    item.prior = NULL;  mono.constraint = mono.constraint;
    conv.crit = 1e-4; maxitr = maxitr; bound = 1e-4
    Ki <- rowSums(Q)
    item.prior <- llply(2^Ki,function(x)matrix(2,nrow = x,ncol = 2))
    for(i in seq_len(I)){
      item.prior[[i]][1,1] <- item.prior[[i]][2^Ki[i],2] <- 1.5
      item.prior[[i]][2^Ki[i],1] <- item.prior[[i]][1,2] <- 2.5

    }
    Ki <- rowSums(Q > 0)
    Li <- 2^Ki
    parloc <- LC2LG(Q=Q)
    prior <- rep(1/L, L)
    logprior <- log(prior)

    item.parm <- list()
    for (i in seq_len(I)) {
      if (Ki[i] == 1) {item.parm[[i]] <- c(0.2, 0.8)}
      else if (Ki[i] == 2) {item.parm[[i]] <- c(0.2, 0.5, 0.5, 0.8)}
      else if (Ki[i] == 3) {item.parm[[i]] <- c(0.2, 0.4, 0.4, 0.4, 0.6, 0.6, 0.6, 0.8)}
      else if (Ki[i] > 3){item.parm[[i]] <- c(0.2, rep(0.5,2^Ki[i]-2), 0.8)}
    }
    item.parm <- l2m(item.parm)

    ConstrMatrix <- vector("list", I)
    ConstrPairs <- vector("list", I)
    for(i in seq_len(I)) {
      ConstrPairs[[i]] <- partial_order2(Ki[i])
      nctj <- nrow(ConstrPairs[[i]])
      tmp <- matrix(0,nctj,2^Ki[i])
      tmp[matrix(c(seq_len(nctj),ConstrPairs[[i]][,1]),ncol = 2)] <- 1
      tmp[matrix(c(seq_len(nctj),ConstrPairs[[i]][,2]),ncol = 2)] <- -1
      ConstrMatrix[[i]] <- tmp
    }

    itr <- 0L
    parm0 <- c(item.parm)
    success <- TRUE
    while(itr < maxitr){
      estep <- LikNR(as.matrix(item.parm), as.matrix(Y), as.matrix(logprior), rep(1,N), as.matrix(parloc), rep(1,N), TRUE)
      Rg = estep$Rg; Ng = estep$Ng

      for(i in 1:I){
        Nj=Ng[i, 1:2^Ki[i]]; Rj=Rg[i, 1:2^Ki[i]]
        r1 <- c(item.prior[[i]][,1]) - 1
        r2 <- c(item.prior[[i]][,2]) - 1
        n <- r1 + r2

        if(any((Nj + n)==0)){
          warning(paste("Nj contains 0 for item", i),call. = FALSE)
          cat("\nFor item ", i, "\n")
          cat(data.frame(Rj=Rj,Nj=Nj,Pj=Rj/Nj))
          return(list(success=FALSE))
        }
        Pj <- (Rj + r1)/(Nj + n)
        Pj[Pj <= bound] <- bound
        Pj[Pj >= 1 - bound] <- 1 - bound

        if(mono.constraint&&any(c(ConstrMatrix[[i]]%*%Pj)<0)){

          obj <- function(x0){-1*sum(Rj*log(x0)+(Nj-Rj)*log(1-x0))-sum(r1*log(x0)+r2*log(1-x0))}
          dev <- function(x0){-1*(Rj + r1 - (Nj + n) * x0)/(x0-x0^2)}
          ineq <- function(x0){c(ConstrMatrix[[i]] %*% x0)}
          ineq.jac <- function(x0){ConstrMatrix[[i]]}
          x00 <- item.parm[i, 1:(2^Ki[i])]

          x00[x00<bound] <- bound
          x00[x00>1-bound] <- 1-bound

          optims <- try(slsqp(x0 = x00,fn=obj, gr = dev, hin=ineq,
                                      hinjac = ineq.jac, lower = rep(bound,length(x00)),
                                      upper = rep(1-bound,length(x00))),silent = TRUE)

          if(inherits(optims,"try-error")){
            warning(paste("Optimization failed for item", i),call. = FALSE)
            if(verbose==1)
              cat(optims)
            return(list(success=FALSE))
          }

          item.parm[i, 1:(2^Ki[i])] <- optims$par
        }else{
          item.parm[i, 1:(2^Ki[i])] <- Pj
        }

      }
      prior <- c(exp(estep$logprior))
      prior <- prior/sum(prior)
      logprior <- log(prior)
      parm1 <- c(item.parm)
      maxchg = max(abs(parm1-parm0),na.rm = TRUE)

      parm0 <- parm1
      itr <- itr + 1

      if(is.infinite(estep$LL))
        stop("-2LL is not finite.",call. = FALSE)
      if(verbose==1L) {
        cat('\rIter =',itr,' Max. abs. change =',formatC(maxchg,digits = 5, format = "f"),
            ' Deviance  =',formatC(-2 * estep$LL,digits = 3, format = "f"),'                                                                                 ')
      }else if (verbose==2L) {
        cat('Iter =',itr,' Max. abs. change =',formatC(maxchg,digits = 5, format = "f"),
            ' Deviance  =',formatC(-2 * estep$LL,digits = 3, format = "f"),'                                                                                \n')
      }

      if(maxchg < conv.crit) break
    }
    estep <- LikNR(as.matrix(item.parm), as.matrix(Y), as.matrix(logprior),
                   rep(1,N), as.matrix(parloc), rep(1,N), FALSE)

    CDM.obj <- list(catprob.parm = m2l(item.parm), posterior.prob = exp(estep$logprior),
                    success=success,item.prior=item.prior)
    CDM.obj$Y <- Y
    CDM.obj$Q <- Q
    CDM.obj$testfit$Deviance <- -2 * estep$LL
    CDM.obj$testfit$npar <- sum(2^rowSums(Q)) - I + L - 1
    CDM.obj$testfit$AIC <- CDM.obj$testfit$Deviance + 2*CDM.obj$testfit$npar
    CDM.obj$testfit$BIC <- CDM.obj$testfit$Deviance + log(N)*CDM.obj$testfit$npar
    P.alpha.Xi <- exp(estep$logpost)
    P.alpha <- colSums(P.alpha.Xi) / sum(P.alpha.Xi)
    alpha.P <- t(apply(P.alpha.Xi, 1, function(x)(return(colSums(x*KS)))))
    alpha <- KS[apply(P.alpha.Xi, 1, which.max), ]
    if(verbose != 0)
      cat("\n")
  }

  res <- list(analysis.obj=CDM.obj, alpha=alpha, P.alpha.Xi=P.alpha.Xi,
              alpha.P=alpha.P, P.alpha=P.alpha,
              model.fit=data.frame(Deviance=CDM.obj$testfit$Deviance,
                                   npar=CDM.obj$testfit$npar,
                                   AIC=CDM.obj$testfit$AIC,
                                   BIC=CDM.obj$testfit$BIC))
  class(res) <- "CDM.obj"

  return(res)
}
