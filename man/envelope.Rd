\name{envelope}
\alias{envelope}
\title{Building of normal probability plots with simulated envelope of the deviance-type residuals.}
\description{
\bold{envelope} is used to calculate and display normal probability plots with simulated envelope of the deviance-type residuals.}
\usage{
envelope(object, reps, conf)
}
\arguments{
  \item{object}{an object of the class \code{ssym}. This object is obtained from the application of \code{ssym.l()} or \code{ssym.nl()}.}
  \item{reps}{a positive integer representing the number of iterations in which the simulated envelopes are based. Default is \code{reps=25}.}
  \item{conf}{value within the interval (0,1) that represents the confidence level of the simulated envelopes. Default is \code{conf=0.95}.}
}

\details{
Objects of the class \code{ssym} obtained from the application of \code{ssym.l2()} are not supported. The smoothing parameters are assumed to be known.
}

\references{Atkinson, A. C. (1985) Plots, transformations and regression: an introduction to graphical methods of diagnostic regression analysis. Oxford Science Publications, Oxford.}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\examples{
###################################################################################
################# Blood flow Data - a log-power-exponential model #################
###################################################################################
#data("la", package="gamlss.nl")
#fit <- ssym.nl(log(PET60) ~ log(bflow) + log(1+b1*exp(-b2/bflow)) | bflow,
#       data=la, start=c(b1=-0.6,b2=98), family="Powerexp", xi=-0.45)
#summary(fit)
#
################## Simulated envelopes ##################
#envelope(fit,reps=50,conf=0.99)
 
}

