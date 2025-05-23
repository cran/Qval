# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

P_GDINA <- function(Qi, P_est, pattern, P_alpha) {
    .Call(`_Qval_P_GDINA`, Qi, P_est, pattern, P_alpha)
}

beta_Ni_ri <- function(pattern, AMP, Y) {
    .Call(`_Qval_beta_Ni_ri`, pattern, AMP, Y)
}

calculatePEst <- function(Yi, P_alpha_Xi) {
    .Call(`_Qval_calculatePEst`, Yi, P_alpha_Xi)
}

get_Pattern <- function(s, alpha) {
    .Call(`_Qval_get_Pattern`, s, alpha)
}

log_likelihood_i <- function(Yi, P_Xj_alpha, P_alpha_Xi) {
    .Call(`_Qval_log_likelihood_i`, Yi, P_Xj_alpha, P_alpha_Xi)
}

