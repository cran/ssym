\name{envelope}
\alias{envelope}
\title{Building of normal probability plots with simulated envelope of the deviance-type residuals.}
\description{
\bold{envelope} is used to calculate and display normal probability plots with simulated envelope of the deviance-type residuals.}
\usage{
envelope(object, reps, conf, xlab.mu, ylab.mu, main.mu, xlab.phi, ylab.phi, main.phi)
}
\arguments{
  \item{object}{an object of the class \code{ssym}. This object returned from the call to \code{ssym.l()} or \code{ssym.nl()}.}
  \item{reps}{a positive integer representing the number of iterations in which the simulated envelopes are based. Default is \code{reps=25}.}
  \item{conf}{value within the interval (0,1) that represents the confidence level of the simulated envelopes. Default is \code{conf=0.95}.}
  \item{xlab.mu}{character. An optional label for the \emph{x} axis for the graph of the deviance-type residuals for the median submodel.}
  \item{ylab.mu}{character. An optional label for the \emph{y} axis for the graph of the deviance-type residuals for the median submodel.}
  \item{main.mu}{character. An optional overall title for the plot for the graph of the deviance-type residuals for the median submodel.}
  \item{xlab.phi}{character. An optional label for the \emph{x} axis for the graph of the deviance-type residuals for the skewness submodel.}
  \item{ylab.phi}{character. An optional label for the \emph{y} axis for the graph of the deviance-type residuals for the skewness submodel.}
  \item{main.phi}{character. An optional overall title for the plot for the graph of the deviance-type residuals for the skewness submodel.}
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

