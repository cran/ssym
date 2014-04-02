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
Date: \tab 2014-03-30\cr
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
library(NISTnls)
Chwirut <- Chwirut1[order(Chwirut1$x,Chwirut1$y),]
attach(Chwirut)

loc.f <- function(vP){
	          beta <- vP[1:3]
	         -beta[1]*x - log(beta[2] + beta[3]*x)	  
}
start <- c(b1=0.15, b2=0.005, b3=0.012)
fit <- ssym.nl(log(y),loc.f,start=start,ncs=x,family='Hyperbolic',xi=1,
               local.influence=TRUE)
summary(fit)
}
