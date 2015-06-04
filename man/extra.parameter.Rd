\name{extra.parameter}
\alias{extra.parameter}
\title{Tool that supports the estimation of the extra parameter.}
\description{
\bold{extra.parameter} is used to plot a graph of the behaviour of the overall goodness-of-fit statistic and \eqn{-2\textsf{L}(\hat{\theta})}
versus the extra parameter \eqn{\zeta} in the interval/region defined by the arguments \code{lower} and \code{upper}.
These graphs may be used to choosing the extra parameter value.}
\usage{
extra.parameter(object, lower, upper)
}
\arguments{
  \item{object}{an object of the class \code{ssym}. This object is obtained from the application of \code{ssym.l()}, \code{ssym.nl()} or \code{ssym.l2()}}
  \item{lower}{lower limit(s) of the interest interval/region for the extra parameter.}
  \item{upper}{upper limit(s) of the interest interval/region for the extra parameter.}
}

\references{Vanegas, L.H. and Paula, G.A. (2015b) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics (to appear)}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
###################################################################################
############### Textures of snacks Data - a log-Student-t model  #################
###################################################################################
#data("Snacks", package="ssym")
#fit <- extra.parameter(log(texture) ~ type + ncs(week) | type, data=Snacks,
#        family='Student', xi=10)
#summary(fit)
#
############################ Extra parameter ###########################
#extra.parameter(fit,5,50)

###################################################################################
################## Biaxial Fatigue Data - a Birnbaum-Saunders model   #############
###################################################################################
#data("Biaxial", package="ssym")
#fit <- ssym.nl(log(Life) ~ b1*Work^b2, start=c(b1=16, b2=-0.25),
#                data=Biaxial, family='Sinh-normal', xi=1.54)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,1.3,1.8)
}

