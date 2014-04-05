\name{lambda.hat}
\alias{lambda.hat}
\title{Choosing the smoothing parameter}
\description{
\bold{lambda.hat} is used to choose the smoothing parameter based on the cross-validation score.}
\usage{
lambda.hat(response,x,xx,lambda,plot)
}
\arguments{
  \item{response}{the response variable.}
  \item{x}{knots of the natural cubic spline, it is a set (ordered ascending) formed by the different values of \emph{t} in the sample.}
  \item{xx}{values of the explanatory variable \emph{t}.}  
  \item{lambda}{a numeric value of starting estimate for the smoothing parameter.}
  \item{plot}{if \emph{plot=TRUE}, it shows a graph with the behavior of the cross-validation score.}  
}
\value{
\item{lh}{a list with the chosen smoothing parameter and value of the cross-validation score evaluated in it.}
}
\references{Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}
\examples{
par(mfrow=c(1,2))
m1 <- "Estimated Natural Cubic Spline"
n <- 300
t <- sort(round(runif(n),digits=2))
y <- cos(4*pi*t) + rnorm(n)
ss <- splinek(as.numeric(levels(factor(t))),t)
N <- ss$N
M <- ss$K
lambda <- lambda.hat(y,as.numeric(levels(factor(t))),t,1,plot=TRUE)
lambda <- lambda$lambda_hat
h <- solve(t(N)\%*\%N + lambda*M)\%*\%t(N)\%*\%y
gam <- solve(ss$R)\%*\%t(ss$Q)\%*\%h
sa <- ncs.graph(as.numeric(levels(factor(t))),h,gam,1000)

plot(t,y,xlim=range(t),ylim=range(y),cex=0.3,lwd=3,xlab="",ylab="",main=m1)
par(new=TRUE)
plot(t,cos(4*pi*t),xlim=range(t),ylim=range(y),type="l",xlab="",ylab="",col="red")
par(new=TRUE)
plot(sa[,1],sa[,2],xlim=range(t),ylim=range(y),type="l",xlab="t",ylab="y",col="blue")
legend(min(t),max(y),col=c("red","blue"),bty="n",lty=1,legend=c("True function","NCS"))
}
