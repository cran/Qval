% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S3Extract.R
\name{extract}
\alias{extract}
\alias{extract.CDM}
\alias{extract.validation}
\alias{extract.sim.data}
\title{Extract Elements from Various Objects}
\usage{
extract(object, what, ...)

\method{extract}{CDM}(object, what, ...)

\method{extract}{validation}(object, what, ...)

\method{extract}{sim.data}(object, what, ...)
}
\arguments{
\item{object}{An object of class \code{\link[Qval]{CDM}}, \code{\link[Qval]{validation}}, or \code{\link[Qval]{sim.data}}.}

\item{what}{A character string naming the component to extract.}

\item{...}{Additional arguments (unused).}
}
\value{
The requested component, whose type depends on \code{what} and the class of \code{object}.
}
\description{
A unified extractor for retrieving components from objects generated by the Qval package.
}
\details{
This generic function and its methods allow users to extract internal elements from
objects of class \code{\link[Qval]{CDM}}, \code{\link[Qval]{validation}}, or \code{\link[Qval]{sim.data}} without directly
accessing object slots. Supported components vary by class:

\describe{
  \item{\emph{\code{\link[Qval]{CDM}}}}{\describe{
    \item{analysis.obj}{The underlying GDINA or BM fit object.}
    \item{alpha}{EAP-based attribute estimates for each individual.}
    \item{P.alpha.Xi}{Posterior probabilities of attribute patterns.}
    \item{alpha.P}{Marginal mastery probability matrix.}
    \item{P.alpha}{Prior attribute weights from final iteration.}
    \item{Deviance}{Model deviance (-2 * marginal log-likelihood).}
    \item{npar}{Number of estimated parameters.}
    \item{AIC}{Akaike Information Criterion.}
    \item{BIC}{Bayesian Information Criterion.}
    \item{call}{Original function call.}
  }}
  \item{\emph{\code{\link[Qval]{validation}}}}{\describe{
    \item{Q.orig}{Original Q-matrix.}
    \item{Q.sug}{Suggested Q-matrix.}
    \item{time.cost}{time cost for validation.}
    \item{process}{Iteration-by-item modification log.}
    \item{iter}{Total iterations performed.}
    \item{priority}{Attribute priority matrix (PAA method only).}
    \item{Hull.fit}{Data for Hull plot (Hull method only).}
    \item{call}{Original function call.}
  }}
  \item{\emph{\code{\link[Qval]{sim.data}}}}{\describe{
    \item{dat}{Simulated response matrix (N × I).}
    \item{Q}{Simulated Q-matrix.}
    \item{attribute}{True attribute pattern matrix (N × K).}
    \item{catprob.parm}{Category success probabilities list.}
    \item{delta.parm}{Delta parameter list.}
    \item{higher.order.parm }{Higher-order parameter list.}
    \item{mvnorm.parm}{Multivariate normal distribution parameters.}
    \item{LCprob.parm}{Item/category success probabilities matrix.}
    \item{call}{Original function call.}
  }}
}
}
\section{Methods (by class)}{
\itemize{
\item \code{extract(CDM)}: Extract fields from a CDM object

\item \code{extract(validation)}: Extract fields from a validation object

\item \code{extract(sim.data)}: Extract fields from a sim.data object

}}
\examples{
library(Qval)
set.seed(123)

\dontrun{
################################################################
# Example 1: sim.data extraction                               #
################################################################
Q <- sim.Q(3, 10)
data.obj <- sim.data(Q, N = 200)
extract(data.obj, "dat")


################################################################
# Example 2: CDM extraction                                    #
################################################################
CDM.obj <- CDM(data.obj$dat, Q)
extract(CDM.obj, "alpha")
extract(CDM.obj, "AIC")


################################################################
# Example 3: validation extraction                             #
################################################################
validation.obj <- validation(data.obj$dat, Q, CDM.obj)
extract(validation.obj, "Q.sug")
}

}
