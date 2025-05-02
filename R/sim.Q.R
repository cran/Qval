#' Generate a Random Q-matrix
#'
#' @description
#' generate a \eqn{I} × \eqn{K} Q-matrix randomly, which consisted of one-attribute q-vectors
#' (50%), two-attribute q-vectors (25%), and three-attribute q-vectors (25%). 
#' This function ensures that the generated Q-matrix contains at least two identity matrices as a priority.
#' Therefore, this function must also satisfy the condition that the number of items (\eqn{I})  
#' must be at least twice the number of attributes (\eqn{K}).  
#'
#' @param I The number of items.
#' @param K The number of attributes of the Q-matrix.
#'
#' @return An object of class \code{matrix}.
#'
#' @author Haijiang Qin <Haijiang133@outlook.com>
#'
#' @references
#' Najera, P., Sorrel, M. A., de la Torre, J., & Abad, F. J. (2021). Balancing fit and parsimony to improve Q-matrix validation. Br J Math Stat Psychol, 74 Suppl 1, 110-130. DOI: 10.1111/bmsp.12228.
#'
#' @examples
#' library(Qval)
#'
#' set.seed(123)
#'
#' Q <- sim.Q(5, 10)
#' print(Q)
#'
#' @export
#' @importFrom GDINA attributepattern
#'

sim.Q <- function(K, I){

  if(K < 2)
    stop("The number of attributes (K) must be more than 1.\n")
  if(I < 2*K)
    stop("The number of items (I) must be twice or more than the number of attributes (K).\n")

  alpha <- attributepattern(K)
  Q.r <- rbind(alpha[2:(K+1), ], alpha[2:(K+1), ])

  I.1 <- ifelse(ceiling(I*0.5) > 2*K, ceiling(I*0.5), 2*K)
  I.2 <- ceiling((I - I.1) / 2)
  I.3 <- I - I.1 - I.2

  Q.1 <- rbind(Q.r, alpha[sample(2:(K+1),
                                 ifelse(I.1 - 2*K > 0, I.1 - 2*K, 0), replace = TRUE), ])
  Q.2 <- alpha[sample((1+choose(K, 1)+1):(1+choose(K, 1)+choose(K, 2)),
                      ifelse(I.2 > 0, I.2, 0), replace = TRUE), ]
  if(K > 3)
    Q.3 <- alpha[sample((1+choose(K, 1)+choose(K, 2)+1):(1+choose(K, 1)+choose(K, 2)+choose(K, 3)),
                        ifelse(I.3 > 0, I.3, 0), replace = TRUE), ]
  if(K == 3)
    Q.3 <- matrix(1, I.3, K)
  if(K == 2)
    Q.3 <- alpha[sample((1+choose(K, 1)+1):(1+choose(K, 1)+choose(K, 2)),
                        ifelse(I.3 > 0, I.3, 0), replace = TRUE), ]

  Q <- rbind(Q.1, Q.2, Q.3)
  Q <- Q[sample(1:I, I), ]
  rownames(Q) <- paste0("item ", 1:I)

  return(Q)
}
