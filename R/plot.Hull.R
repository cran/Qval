
#' 
#' Hull Plot
#' 
#' @description
#' This function can provide the Hull plot. The points suggested by the Hull method are marked in red.
#' 
#' @param x A \code{list} containing all the information needed to plot the Hull plot. 
#'          It can be gotten from the \code{validation} object when \code{method} = \code{"Hull"}.
#' @param i A numeric, which represents the item you want to plot Hull curve. 
#' @param ... Additional arguments to be passed to the plotting function.
#' 
#' 
#' @return None. This function is used for side effects (plotting).
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
#' ############### ESA ###############
#' Hull.obj <- validation(data$dat, MQ, CDM.obj, method = "Hull", search.method = "ESA") 
#' Hull.fit <- Hull.obj$Hull.fit
#' 
#' ## plot Hull curve for item 5
#' plot(Hull.fit, 5)
#' 
#' ############### PAA ###############
#' Hull.obj <- validation(data$dat, MQ, CDM.obj, method = "Hull", search.method = "PAA") 
#' Hull.fit <- Hull.obj$Hull.fit
#' 
#' ## plot Hull curve for item 5
#' plot(Hull.fit, 5)
#' }
#' 
#' 
#' 
#' @export
#' @importFrom graphics plot points text
#' 
plot.Hull <- function(x, i, ...){
  
  number.of.parameters <- x[[i]]$number.of.parameters
  fit.index <- x[[i]]$fit.index
  posi <- x[[i]]$posi
  pattern.criterion <- x[[i]]$pattern.criterion
  pattern <- x[[i]]$pattern
  criter <- x[[i]]$criter
  sug <-  which(fit.index == x[[i]]$sug)
  
  plot(number.of.parameters[posi], fit.index[posi], 
       type = "o", pch = 19, 
       main = "Hull plot", 
       xlab = "Number of Parameters", ylab = criter)
  points(number.of.parameters[-posi], fit.index[-posi], pch = 1)
  points(number.of.parameters[sug], fit.index[sug], col = "red", pch = 19, cex=1)
  
  # Create labels
  labels.pattern.criterion <- c(1, pattern.criterion)
  labels <- rep("[", length(labels.pattern.criterion))
  for (k in 1:ncol(pattern))
    labels <- paste0(labels, pattern[labels.pattern.criterion, k])
  labels <- paste0(labels, rep("]", length(labels.pattern.criterion)))
  
  # Automatically adjust label positions
  for (i in 1:length(number.of.parameters)) {
    # Default pos = 3 (top)
    pos <- 3 
    # Check if close to the top or bottom of the plot
    if (fit.index[i] > max(fit.index) * 0.9) {
      pos <- 1  # If close to the top, place the label below
    } else if (fit.index[i] < min(fit.index) * 1.1) {
      pos <- 3  # If close to the bottom, place the label above
    }
    # Check if close to the left or right edges of the plot
    if (number.of.parameters[i] > max(number.of.parameters) * 0.9) {
      pos <- 2  # If close to the right, place the label to the left
    } else if (number.of.parameters[i] < min(number.of.parameters) * 1.1) {
      pos <- 4  # If close to the left, place the label to the right
    }
    
    text(number.of.parameters[i], fit.index[i], labels = labels[i], pos = pos, cex = 0.8, col = "black")
  }
}
