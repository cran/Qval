.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}

.onLoad <- function(libname, pkgname) {
  if (!requireNamespace("GDINA", quietly = TRUE)) {
    stop("GDINA package is required but not installed. Please install it first.", call. = FALSE)
  }

  GDINA_env <- asNamespace("GDINA")

  l2m <- get("l2m", envir = GDINA_env)
  partial_order2 <- get("partial_order2", envir = GDINA_env)
  LikNR <- get("LikNR", envir = GDINA_env)
  m2l <- get("m2l", envir = GDINA_env)

  assign("l2m", l2m, envir = asNamespace(pkgname))
  assign("partial_order2", partial_order2, envir = asNamespace(pkgname))
  assign("LikNR", LikNR, envir = asNamespace(pkgname))
  assign("m2l", m2l, envir = asNamespace(pkgname))

  # message("Successfully imported unexported functions from GDINA.")
}

StartWelcomeMessage <- function(){

  paste("Qval R Package ",
        "(version ", utils::packageDescription("Qval")$Version,
        "; ",utils::packageDescription("Qval")$Date, ")\n",
        sep="")
}
