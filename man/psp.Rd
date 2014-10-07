\name{psp}
\alias{psp}
\title{Tool to build incidence and penalty matrices in order to approximate smooth functions by natural cubic splines.}
\description{
\bold{psp} is used to build the incidence matrix, the penalty matrix as well as other interest matrices to approximate functions by natural cubic splines.}
\usage{
psp(xx, lambda, nknots, diff)
}
\arguments{
  \item{xx}{values of the explanatory variable.}
  \item{lambda}{an optional positive value, which corresponds to the smoothing parameter.}
  \item{nknots}{an optional integer specifying the number of internal knots. Defaul is \emph{[n^(1/3)]}}.
  \item{diff}{an optional integer specifying the order of the difference penalty term. Default is 2.}  
}
\value{
\item{xx}{a matrix having the following attributes: B-spline basis matrix, penalty matrix, and smoothing parameter (if it has been specified).}
}
\references{Eilers P.H.C. and Marx B.D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science. 11, 89-121.}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
n <- 300
t <- sort(round(runif(n),digits=2))

t2 <- psp(t, diff=3)
N <- attr(t2, "N") ## B-spline basis matrix Matrix
M <- attr(t2, "K") ## Penalty Matrix
}
