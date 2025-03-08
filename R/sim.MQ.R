#' Simulate mis-specifications in the Q-matrix
#'
#' @description
#' simulate certen \code{rate} mis-specifications in the Q-matrix.
#'
#' @param Q The Q-matrix (\code{\link[Qval]{sim.Q}}) that need to simulate mis-specifications.
#' @param rate The ratio of mis-specifications in the Q-matrix.
#' @param verbose Logical indicating to print information or not. Default is \code{TRUE}
#'
#' @return An object of class \code{matrix}.
#'
#' @author Haijiang Qin <Haijiang133@outlook.com>
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' Q <- sim.Q(5, 10)
#' print(Q)
#'
#' MQ <- sim.MQ(Q, 0.1)
#' print(MQ)
#'
#' @export
#' @importFrom GDINA attributepattern
#' @importFrom stats runif
#' 
#'

sim.MQ <- function(Q, rate, verbose = TRUE){
  Q <- as.matrix(Q)
  sum <- floor(ncol(Q) * nrow(Q) * rate)
  sum.1 <- round(runif(1, 0, sum))
  sum.1 <- ifelse(sum.1 <= round(sum(Q)*0.8), sum.1, round(sum(Q)*0.8))
  sum.0 <- sum - sum.1
  wrong.Q <- Q
  if(sum.0 > 0){isWrong.0 <- sample(which(Q == 0), sum.0, replace = FALSE)
  }else{isWrong.0 <- c()}
  if(sum.1 > 0){isWrong.1 <- sample(which(Q == 1), sum.1, replace = FALSE)
  }else{isWrong.1 <- c()}
  wrong.Q[c(isWrong.0, isWrong.1)] <- 1 - wrong.Q[c(isWrong.0, isWrong.1)]

  iter <- 1
  while (any(colSums(wrong.Q) == 0) | any(rowSums(wrong.Q) == 0)) {
    iter <- iter + 1
    if(iter %% 10000 == 0){
      sum.1 <- round(runif(1, 0, sum))
      sum.1 <- ifelse(sum.1 <= round(sum(Q)*0.8), sum.1, round(sum(Q)*0.8))
      sum.0 <- sum - sum.1
    }
    wrong.Q <- Q
    if(sum.0 > 0){isWrong.0 <- sample(which(Q == 0), sum.0, replace = FALSE)
    }else{isWrong.0 <- c()}
    if(sum.1 > 0){isWrong.1 <- sample(which(Q == 1), sum.1, replace = FALSE)
    }else{isWrong.1 <- c()}
    wrong.Q[c(isWrong.0, isWrong.1)] <- 1 - wrong.Q[c(isWrong.0, isWrong.1)]
  }
  if(verbose){
    cat("rate of mis-specifications = ", rate, "\n",
        "rate of  over-specifications = ", round(length(which((Q - wrong.Q) == -1))/(ncol(Q) * nrow(Q)), 2), "\n",
        "rate of under-specifications = ", round(length(which((Q - wrong.Q) == 1))/(ncol(Q) * nrow(Q)), 2), "\n")
  }
  return(wrong.Q)
}
