get.Pattern <- function(s, alpha) {
  pattern <- 0
  for(i in 1:nrow(alpha))
    if(all(s == alpha[i, ])) {
      pattern <- i
      break
    }
  return(pattern)
}
