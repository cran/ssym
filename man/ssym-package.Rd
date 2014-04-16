\name{ssym-package}
\alias{ssym-package}
\alias{ssym}
\docType{package}
\title{Fitting Semiparametric Symmetric Regression Models}
\description{
This package allows to fit a semiparametric regression model suitable for data set analysis in which the conditional distribution of the response is symmetric and continuous.}
\details{
\tabular{ll}{
Package: \tab ssym\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-04-15\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula

Maintainer: Luis Hernando Vanegas
}
\references{Vanegas, L.H. and Paula, G.A. (2014) A Semiparametric Approach for Joint Modeling of Median and Skewness. TEST (submitted)

            Vanegas, L.H. and Paula, G.A. (2014) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics (submitted)}

\examples{
data(Snacks)
Snacks2 <- Snacks[order(Snacks$snack,Snacks$week,Snacks$texture),]
attach(Snacks2)

snack <- factor(snack)

fit <- ssym.l(log(texture), formula.mu=~snack, ncs.mu=week, formula.phi=~snack,
              family='Student', xi=14, local.influence=TRUE)
summary(fit)
}
