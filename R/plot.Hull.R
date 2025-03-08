#' 
#' Hull Plot
#' 
#' @description
#' This function can provide the Hull plot. The point suggested by the Hull method is marked in red.
#' 
#' @param x A \code{\link[Qval]{validation}} in which \code{method} = \code{"Hull"}.
#' @param i A numeric, which represents the item you want to plot Hull curve. 
#' @param ... Additional arguments.
#' 
#' 
#' @return None. This function is used for side effects (plotting).
#' 
#' @examples
#' set.seed(123)
#' library(Qval)
#' 
#' ## generate Q-matrix and data
#' K <- 4
#' I <- 20
#' IQ <- list(
#'   P0 = runif(I, 0.2, 0.4),
#'   P1 = runif(I, 0.6, 0.8)
#' )
#' 
#' \donttest{
#' Q <- sim.Q(K, I)
#' data <- sim.data(Q = Q, N = 500, IQ = IQ, model = "GDINA", distribute = "horder")
#' MQ <- sim.MQ(Q, 0.1)
#' 
#' CDM.obj <- CDM(data$dat, MQ)
#' 
#' ############### ESA ###############
#' Hull.obj <- validation(data$dat, MQ, CDM.obj, method = "Hull", search.method = "ESA") 
#' 
#' ## plot Hull curve for item 20
#' plot(Hull.obj, 20)
#' 
#' ############### PAA ###############
#' Hull.obj <- validation(data$dat, MQ, CDM.obj, method = "Hull", search.method = "PAA") 
#' 
#' ## plot Hull curve for item 20
#' plot(Hull.obj, 20)
#' }
#'  
#' 
#' 
#' @export
#' @importFrom graphics plot points text axis
#' 
plot.validation <- function(x, i, ...){
  
  if(is.null(x$Hull.fit)){
    stop("can not plot Hull when method != 'Hull'")
  }
  Hull.fit <- x$Hull.fit
  
  number.of.parameters <- Hull.fit[[i]]$number.of.parameters
  fit.index <- Hull.fit[[i]]$fit.index
  posi <- Hull.fit[[i]]$posi
  pattern.criterion <- Hull.fit[[i]]$pattern.criterion
  pattern <- Hull.fit[[i]]$pattern
  criter <- Hull.fit[[i]]$criter
  sug <-  which(fit.index == Hull.fit[[i]]$sug)
  
  plot(number.of.parameters[posi], fit.index[posi], 
       type = "o", pch = 19, 
       main = "Hull plot", 
       xlab = "Number of Parameters", ylab = criter, 
       xaxt = "n")
  axis(1, at = number.of.parameters, labels = number.of.parameters)
  points(number.of.parameters[-posi], fit.index[-posi], pch = 1)
  points(number.of.parameters[sug], fit.index[sug], col = "red", pch = 19, cex=1)
  
  # Create labels
  labels.pattern.criterion <- c(1, pattern.criterion)
  labels <- rep("[", length(labels.pattern.criterion))
  for (k in 1:ncol(pattern))
    labels <- paste0(labels, pattern[labels.pattern.criterion, k])
  labels <- paste0(labels, rep("]", length(labels.pattern.criterion)))
  
  size = 1.0
  num.points <- length(number.of.parameters)
  text(number.of.parameters[1], fit.index[1], labels = labels[1], pos = 4, cex = size, col = "black")
  text(number.of.parameters[num.points]*0.975, fit.index[num.points], labels = labels[num.points], pos = 1, cex = size, col = "black")
  
  # Automatically adjust label positions
  for (i in (num.points-1):2) {
    # Default pos = 3 (top)
    pos <- 3
    # Check if close to the top or bottom of the plot
    if (fit.index[i] > max(fit.index) * 0.95) {
      pos <- 1  # If close to the top, place the label below
    }
    
    text(number.of.parameters[i], fit.index[i], labels = labels[i], pos = pos, cex = size, col = "black")
  }
}
