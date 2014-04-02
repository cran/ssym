\name{ncs.graph}
\alias{ncs.graph}
\title{Tool for plotting a natural cubic spline}
\description{
\bold{ncs.graph} is used to find the value of the natural cubic spline \emph{f(t)} at any point \emph{t}, and hence to plot a graph of \emph{f(t)} to any desired degree of accuracy.}
\usage{
ncs.graph(x,g,gam,cont)
}
\arguments{
  \item{x}{knots of the natural cubic spline, it is a set (ordered ascending) formed by the different values of \emph{t} in the sample.}
  \item{g}{a vector with the estimates of \emph{f(t)} for all \emph{t} in \emph{x}.}
  \item{gam}{a vector with the second derivative of \emph{f(t)} for all \emph{t} in \emph{x}.}
  \item{cont}{a positive integer, which indicates the degree of accuracy to plot a graph of \emph{f(t)}.}
}
\value{
\item{pg}{a matrix with two columns representing an ordered set of pairs \emph{(t,f(t))} to plot a graph of \emph{f(t)} to desired degree of accuracy.}
}
\references{Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}
\examples{
m1 <- "Estimated NCS for different values of the smoothing parameter (k)"
n <- 300
t <- round(runif(n),digits=2)
y <- cos(4*pi*t) + rnorm(n)
ss <- splinek(as.numeric(levels(factor(t))),t)
N <- ss$N
M <- ss$K
lambda <- c(0.002,0.02,0.2,2,20)
plot(t,y,xlim=range(t),ylim=range(y),cex=0.3,lwd=3,xlab="",ylab="",main=m1)
for(i in 1:length(lambda)){
	h <- solve(t(N)\%*\%N + lambda[i]*M)\%*\%t(N)\%*\%y
	gam <- solve(ss$R)\%*\%t(ss$Q)\%*\%h
	sa <- ncs.graph(as.numeric(levels(factor(t))),h,gam,1000)
	par(new=TRUE)
	plot(sa[,1],sa[,2],xlim=range(t),ylim=range(y),type="l",xlab="t",ylab="",col=i)
}
legend(min(t),max(y),bty="n",col=1:length(lambda),legend=paste("k = ",lambda),lty=1)
}
