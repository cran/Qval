% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta.R
\name{parallel_iter}
\alias{parallel_iter}
\title{A tool for the \eqn{\beta} Method}
\usage{
parallel_iter(
  i,
  Y,
  P.alpha.Xi,
  P.alpha,
  pattern,
  ri,
  Ni,
  Q.pattern.ini,
  model,
  criter,
  search.method,
  P_GDINA,
  Q.beta,
  L,
  K,
  alpha.P,
  get.MLRlasso,
  priority
)
}
\arguments{
\item{i}{Item number that need to be validated.}

\item{Y}{Observed data matrix for validation.}

\item{P.alpha.Xi}{Individual posterior}

\item{P.alpha}{Attribute prior weights.}

\item{pattern}{The attribute mastery pattern matrix.}

\item{ri}{A vector that contains the numbers of examinees in each knowledge state who correctly answered item \eqn{i}.}

\item{Ni}{A vector that contains the total numbers of examinees in each knowledge state.}

\item{Q.pattern.ini}{Initial pattern number for the model.}

\item{model}{Model object used for fitting (e.g., GDINA).}

\item{criter}{Fit criterion ("AIC", "BIC", "CAIC", or "SABIC").}

\item{search.method}{Search method for model selection ("beta", "ESA", "SSA", or "PAA").}

\item{P_GDINA}{Function to calculate probabilities for GDINA model.}

\item{Q.beta}{Q-matrix for validation.}

\item{L}{Number of latent pattern.}

\item{K}{Number of attributes.}

\item{alpha.P}{Individuals' marginal mastery probabilities matrix (Tu et al., 2022)}

\item{get.MLRlasso}{Function for Lasso regression with multiple linear regression.}

\item{priority}{Vector of priorities for PAA method search.}
}
\value{
An object of class \code{validation} is a \code{list} containing the following components:
\item{fit.index.pre}{The previous fit index value after applying the selected search method.}
\item{fit.index.cur}{The current fit index value after applying the selected search method.}
\item{Q.pattern.cur}{The pattern that corresponds to the optimal model configuration for the current iteration.}
\item{priority}{The priority vector used in the PAA method, if applicable.}
}
\description{
This function performs a single iteration of the \eqn{\beta} method for A item's validation. It is designed 
to be used in parallel computing environments to speed up the validation process of the \eqn{\beta} method. 
The function is a utility function for \code{\link[Qval]{validation}}, and it should not be called independently by the user.
}
