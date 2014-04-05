\name{Snacks}
\alias{Snacks}
\docType{data}
\title{Texture of five different forms of Snacks}
\description{
This data set comes from an experiment developed in the School of Public Health (University of Sao Paulo), in which four different forms of light snacks
(named B, C, D, and E) were compared across 20 weeks with a traditional snack (named A). For the light snacks, the hydrogenated vegetable fat (hvf) was
replaced by canola oil under different proportions: B (0% hvf, 22% canola oil), C (17% hvf, 5% canola oil), D (11% hvf, 11% canola oil) and E (5% hvf, 17% canola oil),
whereas A (22% hvf, 0% canola oil). The experiment was conducted so that in each even week a random sample of 15 units of each snack type was analyzed
in a laboratory and various variables were measured.
}
\usage{data(Snacks)}
\format{
  A data frame with 750 observations on the following 3 variables.
  \describe{
    \item{\code{texture}}{texture of snack.}
    \item{\code{snack}}{a factor with levels 1-5 which correspond to A-E forms of snacks.}
    \item{\code{week}}{week in which the snack was analyzed.}
  }
}
\references{
Paula, G.A., de Moura, A.S., Yamaguchi, A.M. (2004) Sensorial stability of snacks with canola oil and hydrogenated vegetable fat. Technical Report. Center of
Applied Statistics, University of Sao Paulo (in Portuguese).}
\examples{
data(Snacks)
attach(Snacks)
boxplot(log(texture)~snack,xlab="Form of Snacks",ylab="Log(texture)")
}
