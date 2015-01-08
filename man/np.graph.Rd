\name{np.graph}
\alias{np.graph}
\title{Tool for plotting natural cubic splines and P-splines}
\description{
\bold{np.graph} displays a graph of a fitted non-parametric effect (either natural cubic spline or P-spline) from a object of class ``ssym".}
\usage{
np.graph(object, which, exp, xlab, ylab, main)
}
\arguments{
  \item{object}{object of class ``ssym" produced by \emph{ssym.l()} or \emph{ssym.nl()}.}
  \item{which}{an integer indicating the interest submodel, e.g., 1 indicates location submodel, and 2 indicates dispersion submodel.}
  \item{exp}{logical. Should the fitted non-parametric effect be plotted in exponential scale? Default is \emph{FALSE}.}
  \item{xlab}{an optional title for the \emph{x} axis.}
  \item{ylab}{an optional title for the \emph{y} axis.}
  \item{main}{an optional overall title for the plot.}
}
\references{
Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.
Eilers P.H.C. and Marx B.D. (1996). Flexible smoothing with B-splines and penalties. Statistical Science. 11, 89-121.
}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
data("Ovocytes", package="ssym")
fit <- ssym.l(fraction ~ type + psp(time) | type + psp(time), data=Ovocytes,
              family='Powerexp', xi=-0.65)

par(mfrow = c(1,2))
np.graph(fit, which=1, xlab="Time", main="Location")
np.graph(fit, which=2, exp=TRUE, xlab="Time", main="Dispersion")
}
