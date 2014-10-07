\name{Erabbits}
\alias{Erabbits}
\docType{data}
\title{Age and Eye Lens Weight for Rabbits in Australia}
\description{
In this data set, the dry weight of the eye lens was measured for 71 free-living wild rabbits of known age. Eye lens weight tends to vary much less with
environmental conditions than does total body weight, and therefore may be a much better indicator of age.
}
\usage{data(Erabbits)}
\format{
  A data frame with 71 observations on the following 2 variables.
  \describe{
    \item{\code{age}}{age of rabbit in days.}
    \item{\code{wlens}}{dry weight of eye lens in milligrams.}
  }
}
\references{
Dudzinski, M.L. and Mykytowycz, R. (1961) The eye lens as an indicator of age in the wild rabbit in Australia. CSIRO Wildlife Research, 6: 156-159.

Ratkowsky, D. A. (1983). Nonlinear Regression Modelling. Marcel Dekker, New York.

Wei, B. C. (1998). Exponential Family Nonlinear Models. Springer, Singapore.}

\examples{
data("Erabbits", package="ssym")
plot(Erabbits$age, Erabbits$wlens, xlim=range(Erabbits$age), ylim=range(Erabbits$wlens),
     type="p", cex=0.3, lwd=3, ylab="Dry weight of eye lens (in milligrams)",
	 xlab="Age of the animal (in days)")
}
