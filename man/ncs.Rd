\name{ncs}
\alias{ncs}
\title{Tool to build the basis matrix and the penalty matrix of natural cubic splines.}
\description{
\bold{ncs} builds the basis matrix and the penalty matrix to approximate a smooth function using a natural cubic spline.}
\usage{
ncs(xx, lambda, nknots, all.knots)
}
\arguments{
  \item{xx}{the explanatory variable.}
  \item{lambda}{an optional positive value that represents the smoothing parameter value.}
  \item{nknots}{an optional positive integer that represents the number of knots of the natural cubic spline.  Default is \eqn{m=[n^{\frac{1}{3}}]+3}.
                The knots are located at the quantiles of order \eqn{0/(m-1),1/(m-1),\ldots,(m-1)/(m-1)} of xx.}
  \item{all.knots}{logical. If \code{TRUE}, the set of knots and the set of different values of \eqn{xx} coincide. Default is \code{FALSE}.}  
}
\value{
\item{xx}{the explanatory variable \eqn{xx} with the following attributes: set of knots, basis matrix, penalty matrix, smoothing parameter
(if it was specified), and other interest matrices.}
}
\references{
Lancaster, P. and Salkauskas, K. (1986) Curve and Surface Fitting: an introduction. Academic Press, London.
Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.
}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
n <- 300
t <- sort(round(runif(n),digits=1))

t2 <- ncs(t,all.knots=TRUE)
N <- attr(t2, "N") ## Basis Matrix
M <- attr(t2, "K") ## Penalty Matrix
knots <- attr(t2, "knots") ## Set of knots
}