\name{gdp}
\alias{gdp}
\docType{data}
\title{Gross Domestic Product per capita}
\description{
This dataset corresponds to the Gross Domestic Product per capita (current US$) of 190 countries during 2010.}
\usage{data(gdp)}
\format{
  A data frame with 540 observations on the following 2 variables.
  \describe{
    \item{\code{Country}}{Country.}
    \item{\code{gdp2010}}{Gross Domestic Product per capita (current US$).}
  }
}
\references{
World Bank's DataBank website (http://databank.worldbank.org/data/).
}
\examples{
data("gdp", package="ssym")
par(mfrow=c(1,2))
hist(gdp$gdp2010, xlim=range(gdp$gdp2010), ylim=c(0,0.00015), prob=TRUE, breaks=55,
     col="light gray",border="dark gray", xlab="GDP per capita 2010", main="Histogram")
plot(ecdf(gdp$gdp2010), xlim=range(gdp$gdp2010), ylim=c(0,1), verticals=TRUE,
     do.points=FALSE, col="dark gray", xlab="GDP per capita 2010",
	 main="Empirical Cumulative Distribution Function")
}
