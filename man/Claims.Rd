\name{Claims}
\alias{Claims}
\docType{data}
\title{Personal Injure Insurance}
\description{
This data set contains information on 540 settled personal injury insurance claims from an Australian insurance company, which is a sample of the original data set. These claims
had legal representation and arose from accidents occurring from January 1998 through to June 1999.
}
\usage{data(Claims)}
\format{
  A data frame with 540 observations on the following 2 variables.
  \describe{
    \item{\code{total}}{amount of paid money by an insurance policy in thousands of Australian dollars.}
    \item{\code{accmonth}}{month of occurrence of the accident coded 103 (January 1998) through to 120 (June 1999).}
    \item{\code{op_time}}{operational time in percentage.}	
  }
}
\references{
de Jong P, Heller GZ. Generalized Linear Models for Insurance Data. Cambridge University Press: Cambridge, England, 2008.

Paula, G.A., Leiva, V., Barros, M. and Liu, S. (2012) Robust statistical modeling using the Birnbaum-Saunders-t distribution applied to insurance distribution, Applied Stochastic Model in Business and Industry, 28:16-34.
}
\examples{
data("Claims", package="ssym")
plot(Claims$op_time, Claims$total, xlim=range(Claims$op_time), ylim=range(Claims$total),
     type="p", cex=0.3, lwd=3, ylab="Amount of paid money",xlab="Operational time",
	 main="Personal Injure Insurance Data")
}
