\name{Ovocytes}
\alias{Ovocytes}
\docType{data}
\title{Fraction of cell volume}
\description{
This data set comes from an experiment comparing the responses of immature and mature goat ovocytes to an hyper-osmotic test. As a compound permeates,
water reenters the cell, and the cell re-expands until the system reaches an osmotic equilibrium. The results are obtained using immature and ovulated
(mature) ovocytes exposed to propanediol, a permeable compound. Then, the cell volume during equilibration is recorded at each time \emph{t}.}

\usage{data(Ovocytes)}
\format{
  A data frame with 161 observations on the following 3 variables.
  \describe{
    \item{\code{type}}{stage of the goat ovocyte: \emph{Mature} or \emph{Immature}.}
    \item{\code{time}}{time since exposition to propanediol.}
    \item{\code{fraction}}{fraction of initial isotonic cell volume at any given time \emph{t} during equilibration.}
  }
}
\references{
Huet, S., Bouvier, A., Gruet, M.A. and Jolivet, E. (1996). Statistical Tools for Nonlinear Regression. Springer, New York.

Le Gal F., Gasqui P., Renard J.P. (1994) Differential Osmotic Behavior of Mammalian Oocytes before and after Maturation: A Quantitative Analysis Using
Goat Oocytes as a Model. Cryobiology, 31: 154-170.

Huet S., Bouvier A., Gruet M.A., Jolivet E. (1996) Statistical Tools for Nonlinear Regression. Springer-Verlag: New York.}
\examples{
data("Ovocytes", package="ssym")
xl <- "Time"
yl <- "Fraction of Cell Volume"
mm <- "Fraction of Cell Volume for Mature and Immature Goat Ovocytes"
rx <- range(Ovocytes$time)
ry <- range(Ovocytes$fraction)
plot(Ovocytes$time[Ovocytes$type=='Mature'], Ovocytes$fraction[Ovocytes$type=='Mature'],
     xlim=rx, ylim=ry, type="p", cex=0.5, lwd=1, ylab="", xlab="")
par(new=TRUE)
plot(Ovocytes$time[Ovocytes$type=='Immature'], Ovocytes$fraction[Ovocytes$type=='Immature'],
     xlim=rx, ylim=ry, type="p", cex=0.5, lwd=2, ylab=yl, xlab=xl, main=mm)
legend(rx[1], ry[2], pt.lwd=c(1,2), bty="n", legend=c("Mature","Immature"), pt.cex=0.5, pch=1)
}