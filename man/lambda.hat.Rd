\name{lambda.hat}
\alias{lambda.hat}
\title{Choosing the smoothing parameter}
\description{
\bold{lambda.hat} is used to choose the smoothing parameter using the cross-validation score.}
\usage{
lambda.hat(response,xx,lambda,type,plot)
}
\arguments{
  \item{response}{the response variable.}
  \item{xx}{values of the explanatory variable.}  
  \item{lambda}{an optional numeric value of starting estimate for the smoothing parameter. Default is 1.}
  \item{type}{an integer indicating the type of nonparametric effect, e.g., 1 indicates natural cubic spline, and 2 indicates P-spline.}  
  \item{plot}{logical. If \emph{plot=TRUE}, it shows a graph with the behavior of the cross-validation score.}  
}
\value{
\item{lh}{a list with the chosen smoothing parameter and value of the cross-validation score evaluated in it.}
}
\references{Green, P.J. and Silverman, B.W. (1994) Nonparametric Regression and Generalized Linear Models, Boca Raton: Chapman and Hall.}
\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}
\examples{
n <- 300
t <- sort(round(runif(n),digits=2))
y <- cos(4*pi*t) + rnorm(n)
datas <- data.frame(t,y)

par(mfrow=c(1,2))
lambda <- lambda.hat(y,t,1,1,plot=TRUE)
lambda <- lambda$lambda_hat
fit <- ssym.l(y ~ ncs(t, lambda=lambda), family="Normal", data=datas)
np.graph(fit, which=1, main="Natural Cubic Spline")

#x11()
par(mfrow=c(1,2))
lambda <- lambda.hat(y,t,1,2,plot=TRUE)
lambda <- lambda$lambda_hat
fit2 <- ssym.l(y ~ psp(t, lambda=lambda), family="Normal", data=datas)
np.graph(fit2, which=1, main="P-spline")

}
