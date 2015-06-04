\name{np.graph}
\alias{np.graph}
\title{Tool to plot natural cubic splines or P-splines.}
\description{
\bold{np.graph} displays a graph of a fitted nonparametric effect, either natural cubic spline or P-spline, from an object of class \code{ssym}.}
\usage{
np.graph(object, which, var, exp, simul, obs, xlab, ylab, xlim, ylim, main)
}
\arguments{
  \item{object}{an object of the class \code{ssym}. This object is obtained from the application of \code{ssym.l()}, \code{ssym.nl()} or \code{ssym.l2()}.}
  \item{which}{an integer indicating the interest submodel. For example, 1 indicates location submodel, and 2 indicates skewness (or relative dispersion) submodel.}
  \item{var}{character. It allows to choosing the nonparametric effect using the name of the associated explanatory variable.}
  \item{exp}{logical. If \code{TRUE}, the fitted nonparametric effect is plotted in exponential scale. Default is \code{FALSE}.}
  \item{simul}{logical. If \code{TRUE}, the fitted nonparametric effect is plotted jointly with their \eqn{95\%} simultaneous confidence intervals. If \code{TRUE},
               the fitted nonparametric effect is plotted jointly with their \eqn{95\%} pointwise confidence intervals. Default is \code{TRUE}.}
  \item{obs}{logical. If \code{TRUE}, the fitted nonparametric effect is plotted jointly with the observed data. Default is \code{FALSE}.}
  \item{xlab}{character. An optional label for the \emph{x} axis.}
  \item{ylab}{character. An optional label for the \emph{y} axis.}
  \item{xlim}{numeric. An optional range of values for the \emph{x} axis.}
  \item{ylim}{numeric. An optional range of values for the \emph{y} axis.}
  \item{main}{character. An optional overall title for the plot.}
}
\references{
Lancaster, P. and Salkauskas, K. (1986) Curve and Surface Fitting: an introduction. Academic Press, London.
Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.
Eilers P.H.C. and Marx B.D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science. 11, 89-121.
}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
#data("Ovocytes", package="ssym")
#fit <- ssym.l(fraction ~ type + psp(time) | type + psp(time), data=Ovocytes,
#              family='Powerexp', xi=-0.55)
#
#par(mfrow = c(1,2))
#np.graph(fit, which=1, xlab="Time", main="Location")
#np.graph(fit, which=2, exp=TRUE, xlab="Time", main="Dispersion")
}
