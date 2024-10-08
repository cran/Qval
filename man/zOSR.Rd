% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Qvalindex.R
\name{zOSR}
\alias{zOSR}
\title{Caculate over-specifcation rate (OSR)}
\usage{
zOSR(Q.true, Q.sug)
}
\arguments{
\item{Q.true}{The true Q-matrix.}

\item{Q.sug}{The Q-matrix that has being validated.}
}
\value{
A numeric (OSR index).
}
\description{
Caculate over-specifcation rate (OSR)
}
\details{
The OSR is defned as:
\deqn{
 OSR = \frac{\sum_{i=1}^{I}\sum_{k=1}^{K}I(q_{ik}^{t} < q_{ik}^{s})}{I × K}
}
where \eqn{q_{ik}^{t}} denotes the \code{k}th attribute of item \code{i} in the true Q-matrix (\code{Q.true}),
\eqn{q_{ik}^{s}} denotes \code{k}th attribute of item \code{i} in the suggested Q-matrix(\code{Q.sug}),
and \eqn{I(\cdot)} is the indicator function.
}
\examples{

library(Qval)

set.seed(123)

example.Q1 <- sim.Q(5, 30)
example.Q2 <- sim.MQ(example.Q1, 0.1)
OSR <- zOSR(example.Q1, example.Q2)
print(OSR)

}
