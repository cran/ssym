\name{splinek}
\alias{splinek}
\title{Tool to approximate functions by natural cubic splines}
\description{
\bold{splinek} is used to build the incidence matrix, the penalty matrix as well as other interest matrices to approximate functions by natural cubic splines.}
\usage{
splinek(xx)
}
\arguments{
  \item{xx}{values of the explanatory variable \emph{t}.}  
}
\value{
\item{ss}{a list with the incidence matrix, the penalty matrix as well as other interest matrices.}
}
\references{Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}
