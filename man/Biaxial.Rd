\name{Biaxial}
\alias{Biaxial}
\docType{data}
\title{Brown and Miller's Biaxial Fatigue}
\description{This data set describes the life of a metal piece in cycles to failure. The response is the number of cycles to failure and the explanatory variable is the work per cycle.}
\usage{data(Biaxial)}
\format{
  A data frame with 46 observations on the following 2 variables.
  \describe{
    \item{\code{Work}}{work per cycle.}
    \item{\code{Life}}{number of cycles to failure.}
  }
}
\references{
J.R. Rieck and J.R. Nedelman (1991) A log-linear model for the Birnbaum-Saunders distribution, Technometrics 33, 51:60.}

\examples{
data("Biaxial", package="ssym")
plot(Biaxial$Work, Biaxial$Life, type="p",
     cex=0.3, lwd=3, ylab="Life", xlab="Work per cycle",
	 main="Brown and Miller's Biaxial Fatigue Data")
}
