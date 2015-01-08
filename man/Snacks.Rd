\name{Snacks}
\alias{Snacks}
\docType{data}
\title{Textures of five different types of snacks}
\description{
This data set comes from an experiment developed in the School of Public Health - University of S\~ao Paulo, in which four different forms of light snacks
(denoted by B, C, D, and E) were compared with a traditional snack (denoted by A) for 20 weeks.
For the light snacks, the hydrogenated vegetable fat (hvf) was
replaced by canola oil using different proportions: B (0\% hvf, 22\% canola oil), C (17\% hvf, 5\% canola oil), D (11\% hvf, 11\% canola oil) and E (5\% hvf, 17\% canola oil);
A (22\% hvf, 0\% canola oil) contained no canola oil.
The experiment was conducted such that a random sample of 15 units of each snack type was analyzed
in a laboratory in each even week to measure various variables.
A total of 75 units was analyzed in each even week; with 750 units being analyzed during the experiment.
}
\usage{data(Snacks)}
\format{
  A data frame with 750 observations on the following 3 variables.
  \describe{
    \item{\code{texture}}{texture of the snack unit.}
    \item{\code{type}}{a factor with levels 1-5 which correspond to A-E types of snacks.}
    \item{\code{week}}{week in which the snack unit was analyzed.}
  }
}
\references{
Paula, G.A., de Moura, A.S., Yamaguchi, A.M. (2004) Sensorial stability of snacks with canola oil and hydrogenated vegetable fat. Technical Report. Center of
Applied Statistics, University of Sao Paulo (in Portuguese).

Paula, G.A. (2013) On diagnostics in double generalized linear models. Computational Statistics and Data Analysis, 68: 44-51.}
\examples{
data("Snacks", package="ssym")
boxplot(log(Snacks$texture) ~ Snacks$type, xlab="Type of Snack", ylab="Log(texture)")
}
