\name{ssym-package}
\alias{ssym-package}
\alias{ssym}
\docType{package}
\title{Fitting Semiparametric Log-symmetric Regression Models}
\description{
This package allows to fit a semi-parametric regression model suitable for 
analysis of data sets in which the response variable is continuous, strictly
positive, asymmetric and possibly, censored.}
\details{
\tabular{ll}{
Package: \tab ssym\cr
Type: \tab Package\cr
Version: \tab 1.5.6\cr
Date: \tab 2016-07-21\cr
License: \tab GPL-2 | GPL-3\cr
}
}
\author{
Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula

Maintainer: Luis Hernando Vanegas
}
\references{Vanegas, L.H. and Paula, G.A. (2015) A semiparametric approach for joint modeling of median and skewness. TEST 24, 110-135.

			Vanegas, L.H. and Paula, G.A. (2016) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics 30, 196-220.

			Vanegas, L.H. and Paula, G.A. (2016) An extension of log-symmetric regression models: R codes and applications. Journal of Statistical Computation and Simulation 86, 1709-1735.}

\examples{
data("Snacks", package="ssym")                                    
fit <- ssym.l(log(texture) ~ type + ncs(week) | type, data=Snacks,
               family='Student', xi=15)     
summary(fit)                                                      
}
