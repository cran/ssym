\name{ncs}
\alias{ncs}
\title{Tool to build incidence and penalty matrices in order to approximate smooth functions by natural cubic splines.}
\description{
\bold{ncs} is used to build the incidence matrix, the penalty matrix as well as other interest matrices to approximate functions by natural cubic splines.}
\usage{
ncs(xx, lambda)
}
\arguments{
  \item{xx}{values of the explanatory variable.}
  \item{lambda}{an optional positive value, which corresponds to the smoothing parameter.}
}
\value{
\item{xx}{a matrix having the following attributes: incidence matrix, penalty matrix, smoothing parameter (if it has been specified), and other interest matrices.}
}
\references{Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
n <- 300
t <- sort(round(runif(n),digits=1))

t2 <- ncs(t)
N <- attr(t2, "N") ## Incidence Matrix
M <- attr(t2, "K") ## Penalty Matrix

}