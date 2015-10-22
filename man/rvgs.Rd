\name{rvgs}
\alias{rvgs}
\title{Random generation for some symmetric continuous distributions.}
\description{
\bold{rvgs} is used to random generation from some standard symmetric continuous distributions.}
\usage{
rvgs(n, family, xi)
}
\arguments{
  \item{n}{number of observations.}
  \item{family}{Supported families include \emph{Normal}, \emph{Student}, \emph{Contnormal},
               \emph{Powerexp}, \emph{Hyperbolic}, \emph{Slash}, \emph{Sinh-normal} and \emph{Sinh-t}, which correspond to normal, Student-t, contaminated normal,
			    power exponential, symmetric hyperbolic, slash, sinh-normal and sinh-t distributions, respectively.}
  \item{xi}{a numeric value or numeric vector that represents the extra parameter value of the specified distribution.}
}

\value{
\item{x}{a vector of \eqn{n} observations.}
}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}
\examples{
m1 <- "Standard Sinh-t distributions"
n <- 1000000
xi <- c(10,6,4)
plot(density(rvgs(n,"Sinh-t",xi=c(25,10))), xlim=c(-4.5,4.5), ylim=c(0,0.3), xlab="",
     ylab="", col=1, main="")
par(new=TRUE)
plot(density(rvgs(n,"Sinh-t",xi=c(25,6))), xlim=c(-4.5,4.5), ylim=c(0,0.3), xlab="",
     ylab="", col=2, main="")
par(new=TRUE)
plot(density(rvgs(n,"Sinh-t",xi=c(25,4))), xlim=c(-4.5,4.5), ylim=c(0,0.3), xlab="y",
     ylab="f(y)", main=m1, col=3)
legend(-4, 0.3, bty="n", legend=paste("xi = (",25,",",xi,")"), col=1:4, lty=1)
}
