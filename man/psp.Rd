\name{psp}
\alias{psp}
\title{Tool to build the basis matrix and the penalty matrix of P-splines.}
\description{
\bold{psp} builds the basis matrix and the penalty matrix to approximate a smooth function by using a P-spline.}
\usage{
psp(xx, lambda, b.order, nknots, diff)
}
\arguments{
  \item{xx}{the explanatory variable.}
  \item{lambda}{an optional positive value that represents the smoothing parameter value.}
  \item{b.order}{an optional positive integer that specifies the degree of the B-spline basis matrix. Default is 3.}    
  \item{nknots}{an optional positive integer that represents the number of knots of the P-spline. Default is \eqn{[n^{1/3}]+3}.}
  \item{diff}{an optional positive integer that specifies the order of the difference penalty term. Default is 2.}  
}
\value{
\item{xx}{the explanatory variable \eqn{xx} with the following attributes: set of knots, B-spline basis matrix, penalty matrix
and  smoothing parameter (if it was specified).}
}
\references{Eilers P.H.C. and Marx B.D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science. 11, 89-121.}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
n <- 300
t <- sort(round(runif(n),digits=2))

t2 <- psp(t, diff=3)
N <- attr(t2, "N") ## B-spline basis matrix
M <- attr(t2, "K") ## Penalty Matrix
knots <- attr(t2,"knots") ## Set of knots
}
