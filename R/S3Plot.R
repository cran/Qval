#' @title Plot Methods for Various Qval Objects
#'
#' @description
#' Generate visualizations for objects created by the Qval package. The generic `plot` dispatches to
#' appropriate methods based on object class:
#' \describe{
#'   \item{CDM}{Barplot of attribute-pattern distribution (frequency and proportion).}
#'   \item{sim.data}{Barplot of simulated attribute-pattern distribution (frequency and proportion).}
#'   \item{validation}{Hull plot marking the suggested point in red (method = "Hull").}
#' }
#'
#' @param x An object of class \code{\link[Qval]{CDM}}, \code{\link[Qval]{sim.data}}, or \code{\link[Qval]{validation}}.
#' @param i For \code{validation} objects, the index of the item for which to plot the Hull curve.
#' @param ... Additional arguments (currently unused).
#'
#' @return None. Functions are called for side effects (plotting).
#'
#' @examples
#' set.seed(123)
#' library(Qval)
#'
#' \donttest{
#' K <- 4
#' I <- 20
#' IQ <- list(
#'   P0 = runif(I, 0.2, 0.4),
#'   P1 = runif(I, 0.6, 0.8)
#' )
#' 
#' ################################################################
#' # Example 1: sim.data object                                   #
#' ################################################################
#' Q <- sim.Q(K, I)
#' data.obj <- sim.data(Q = Q, N = 500, IQ = IQ,
#'                      model = "GDINA", distribute = "horder")
#' 
#' plot(data.obj)
#'                      
#' 
#' ################################################################
#' # Example 2: CDM object                                        #
#' ################################################################
#' CDM.obj <- CDM(data.obj$dat, Q, model = "GDINA", 
#'                method = "EM", maxitr = 2000, verbose = 1)
#' plot(CDM.obj)
#'
#'
#' ################################################################
#' # Example 3: validation object (Hull plot)                     #
#' ################################################################
#' MQ <- sim.MQ(Q, 0.1)
#' 
#' CDM.obj <- CDM(data.obj$dat, MQ)
#' 
#' ############### ESA ###############
#' Hull.obj <- validation(data.obj$dat, MQ, CDM.obj, 
#'                        method = "Hull", search.method = "ESA") 
#' 
#' ## plot Hull curve for item 20
#' plot(Hull.obj, 20)
#' 
#' ############### PAA ###############
#' Hull.obj <- validation(data.obj$dat, MQ, CDM.obj, 
#'                        method = "Hull", search.method = "PAA") 
#' 
#' ## plot Hull curve for item 20
#' plot(Hull.obj, 20)
#' }
#'
#' @name plot
NULL

#' @describeIn plot Plot method for CDM objects
#' @importFrom graphics barplot abline text axis mtext par
#' @export
plot.CDM <- function(x, ...) {
  summary.obj <- summary(x)
  patterns <- summary.obj$patterns
  
  patterns.names <- colnames(patterns)
  patterns.freq <- as.numeric(patterns[1, ])
  patterns.prop <- as.numeric(patterns[2, ])
  
  raw_max <- max(patterns.freq) * 1.1
  y_max <- ceiling(raw_max / 10) * 10
  y_ticks <- round(seq(0, y_max, length.out = 10))
  
  par(mar = c(5, 4, 4, 6))
  
  bar_positions <- barplot(patterns.freq, 
                           names.arg = patterns.names,
                           main = paste0("Distribution of Alpha in ", x$analysis.obj$model[1]),
                           xlab = "", 
                           ylab = "Frequency",
                           col = "skyblue",
                           las = 2,
                           ylim = c(0, y_max),
                           yaxt = "n"
  )
  
  abline(h = y_ticks, col = "gray", lty = 2)
  
  text(x = bar_positions, 
       y = patterns.freq + y_max * 0.05, 
       labels = paste0(patterns.freq, "\n(", round(patterns.prop * 100, 1), "%)"), 
       cex = 0.8, col = "black"
  )
  
  axis(side = 2, at = y_ticks, labels = y_ticks, las = 2)
  
  axis(side = 4, 
       at = y_ticks, 
       labels = paste0(round(y_ticks / sum(patterns.freq) * 100, 2), "%"), 
       las = 2, col.axis = "black"
  )
  mtext("Proportion (%)", side = 4, line = 4, col = "black")
}

#' @describeIn plot Plot method for sim.data objects
#' @importFrom graphics barplot abline text axis mtext par
#' @export
plot.sim.data <- function(x, ...) {
  summary.obj <- summary(x)
  patterns <- summary.obj$patterns
  
  patterns.names <- colnames(patterns)
  patterns.freq <- as.numeric(patterns[1, ])
  patterns.prop <- as.numeric(patterns[2, ])
  
  raw_max <- max(patterns.freq) * 1.1
  y_max <- ceiling(raw_max / 10) * 10
  y_ticks <- round(seq(0, y_max, length.out = 10))
  
  par(mar = c(5, 4, 4, 6))
  
  bar_positions <- barplot(patterns.freq, 
                           names.arg = patterns.names,
                           main = "Distribution of Simulated Alpha",
                           xlab = "", 
                           ylab = "Frequency",
                           col = "skyblue",
                           las = 2,
                           ylim = c(0, y_max),
                           yaxt = "n"
  )
  
  abline(h = y_ticks, col = "gray", lty = 2)
  
  text(x = bar_positions, 
       y = patterns.freq + y_max * 0.05, 
       labels = paste0(patterns.freq, "\n(", round(patterns.prop * 100, 1), "%)"), 
       cex = 0.8, col = "black"
  )
  
  axis(side = 2, at = y_ticks, labels = y_ticks, las = 2)
  
  axis(side = 4, 
       at = y_ticks, 
       labels = paste0(round(y_ticks / sum(patterns.freq) * 100, 2), "%"), 
       las = 2, col.axis = "black"
  )
  mtext("Proportion (%)", side = 4, line = 4, col = "black")
}

#' @describeIn plot Hull plot for validation objects
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
  
  par(mar = c(5, 4, 4, 3))
  
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