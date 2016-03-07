\name{myeloma}
\alias{myeloma}
\docType{data}
\title{Survival times for multiple myeloma patients}
\description{The problem is to relate survival times for multiple myeloma patients
to a number of prognostic variables.}
\usage{data("myeloma")}
\format{
  A data frame with 65 observations on the following 7 variables.
  \describe{
    \item{\code{t}}{survival times, in months.}
    \item{\code{event}}{censoring status.}
    \item{\code{x1}}{logarithm of a blood urea nitrogen measurement at diagnosis.}
    \item{\code{x2}}{hemoglobin measurement at diagnosis.}
    \item{\code{x3}}{age at diagnosis.}
    \item{\code{x4}}{sex: 0, male; 1, female.}
    \item{\code{x5}}{serum calcium measurement at diagnosis.}
  }
}
\references{
J.F. Lawless (2002) Statistical Models and Methods for Lifetime Data, Wiley, New York. 
A.P. Li, Z.X. Chen and F.C. Xie (2012) Diagnostic analysis for heterogeneous log-Birnbaum-Saunders regression models, Statistics and Probability Letters 82, 1690:1698.}
\keyword{datasets}
