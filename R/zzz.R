
.onAttach <- function(libname, pkgname){
  packageStartupMessage(StartWelcomeMessage())
}

.onLoad <- function(libname, pkgname) {

  GDINA_env <- asNamespace("GDINA")

  l2m <- get("l2m", envir = GDINA_env)
  partial_order2 <- get("partial_order2", envir = GDINA_env)
  LikNR <- get("LikNR", envir = GDINA_env)
  m2l <- get("m2l", envir = GDINA_env)
  inverse_crossprod <- get("inverse_crossprod", envir = GDINA_env)
  score_pj <- get("score_pj", envir = GDINA_env)
  
  assign("l2m", l2m, envir = asNamespace(pkgname))
  assign("partial_order2", partial_order2, envir = asNamespace(pkgname))
  assign("LikNR", LikNR, envir = asNamespace(pkgname))
  assign("m2l", m2l, envir = asNamespace(pkgname))
  assign("inverse_crossprod", inverse_crossprod, envir = asNamespace(pkgname))
  assign("score_pj", score_pj, envir = asNamespace(pkgname))

}

StartWelcomeMessage <- function(){
  paste("Qval R Package ",
        "(version ", utils::packageDescription("Qval")$Version,
        "; ",utils::packageDescription("Qval")$Date, ")\n",
        sep="")
}

printPackageInfo <- function() {
  packageinfo <- utils::packageDescription("Qval")
  cat(paste("Qval version", packageinfo$Version, "(", packageinfo$Date, ")", sep = ""), "\n")
}
