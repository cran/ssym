\name{ssym.nl}
\alias{ssym.nl}
\title{Fitting Semiparametric Symmetric (or Log-Symmetric) Regression Models}
\description{
\bold{ssym.nl} is used to fit a semiparametric regression model suitable for 
data set analysis in which the conditional distribution of the response (or log-response) is 
symmetric and continuous. In this setup, both location and dispersion
parameters of response (or log-response) variable distribution are explicitly modeled, 
the location using a nonlinear function and the dispersion using a 
semiparametric function, which is aproximated by a natural cubic spline.}
\usage{
ssym.nl(response, mu, start, formula.phi, ncs, start.lambda, lambda, family, xi,
       epsilon, maxiter, subset, local.influence)
}
\arguments{
  \item{response}{the response (or log-response) variable.}
  \item{mu}{an R function with the nonlinear model that describes the location parameter of the response (or log-response) variable distribution.}
  \item{start}{(optional) a named numeric vector of starting estimates.}
  \item{formula.phi}{(optional) a symbolic description of the parametric function to be fitted to the dispersion parameter.}
  \item{ncs}{(optional) an explanatory (continuous) variable to be used in the nonparametric function to be fitted to the dispersion parameter.}
  \item{start.lambda}{(optional) a numeric value of starting estimate for the smoothing parameter.}
  \item{lambda}{(optional) a numerical value for the smoothing parameter indicating that it is provided by the user rather than estimated from the data.}
  \item{family}{a description of the error distribution to be used in the model. Supported families include \emph{Normal}, \emph{Student}, \emph{Powerexp},
                \emph{Hyperbolic}, \emph{Slash}, \emph{Sinh-normal} and \emph{Sinh-t}, which correspond to normal, Student-t, power exponential, symmetric hyperbolic,
				slash, sinh-normal and sinh-t distributions, respectively.}
  \item{xi}{a numeric value or numeric vector that represents the extra parameter of the specified error distribution.}
  \item{epsilon}{(optional) positive convergence tolerance. Default value is 0.0000001.}
  \item{maxiter}{(optional) an positive integer giving the maximal number of iterations for the estimating process. Default value is 500.}
  \item{subset}{(optional) expression indicating individuals or observations to keep (or drop).}
  \item{local.influence}{(optional) logical. If TRUE, local influence measures for the location parameters vector under two perturbation schemes are calculated.}
}
\details{
The iterative estimation process for the parameters of interest is based on the Fisher scoring and backfitting algorithms. Because some distributions such as
Student-t, power exponential, slash and symmetric hyperbolic may be obtained as a scale mixture of normal distributions, the EM algorithm is applied in those cases to
obtain a more efficient iterative process for the parameter estimation. Further, because the Sinh-t distribution can be obtained as a 
scale mixture of Sinh-normal distributions, the EM algorithm is also applied in that case to obtain a more efficient iterative process for the parameter estimation.
The smoothing parameter is chosen using the cross-validation score.
}
\value{
\item{coefs.mu}{a vector of parameter estimates associated with the nonlinear function fitted to the location of the response (or log-response) variable.}
\item{coefs.phi}{a vector of parameter estimates associated with the semiparametric function fitted to the dispersion of the response (or log-response) variable.}
\item{vcov.mu}{approximated variance-covariance matrix asociated with \emph{coefs.mu}.}
\item{se.phi}{approximated standard errors asociated with \emph{coefs.phi}.}
\item{v}{final weights of the iterative process.}
\item{lambda}{estimate for the smoothing parameter.}
\item{gle}{degrees of freedom associated with the nonparametric part of the model.}
\item{deviance.mu}{a vector of \emph{deviances} associated with the location of the response (or log-response) variable.}
\item{deviance.phi}{a vector of \emph{deviances} associated with the dispersion of the response (or log-response) variable.}
\item{mu.fitted}{a vector of fitted values for the location of the response (or log-response) variable.}
\item{phi.fitted}{a vector of fitted values for the dispersion of the response (or log-response) variable.}
\item{lpdf}{a vector of individual contributions to the log-likelihood function.}
\item{pdfz}{a vector of the cumulative distribution function of each individual.}
\item{cw}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the case-weight perturbation scheme)
          asociated with \emph{coefs.mu} and based on the conformal normal curvature.}
\item{pr}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the response perturbation scheme)
          asociated with \emph{coefs.mu} and based on the conformal normal curvature.}
}

\references{Vanegas, L.H. and Paula, G.A. (2014) A Semiparametric Approach for Joint Modeling of Median and Skewness. TEST (submitted)

            Vanegas, L.H. and Paula, G.A. (2014) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics (submitted)}


\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\seealso{
\code{\link{ssym.l}}
}
\examples{
#########################################################################################
########## Ultrasonic Calibration Data - a log-hyperbolic semiparametric model ##########
#########################################################################################
library(NISTnls)
Chwirut <- Chwirut1[order(Chwirut1$x,Chwirut1$y),]
attach(Chwirut)

loc.f <- function(vP){
	          beta <- vP[1:3]
	         -beta[1]*x - log(beta[2] + beta[3]*x)	  
}
start <- c(b1=0.15, b2=0.005, b3=0.012)
fit <- ssym.nl(log(y),loc.f,start=start,ncs=x,family='Hyperbolic',xi=1,
               local.influence=TRUE)
summary(fit)

####################### Plot of the fitted model ##############################

par(mfrow=c(4,2))
xl <- "Distance to the metal"
rx <- range(x)
plot(x,y,xlim=rx,ylim=range(y),type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(x,exp(fit$mu.fitted),xlim=rx,ylim=range(y),type="l",ylab="",xlab=xl,main="Median")

h <- fit$coefs.phi
ss <- splinek(as.numeric(levels(factor(x))),x)
gam <- solve(ss$R)%*%t(ss$Q)%*%h
sa <- ncs.graph(as.numeric(levels(factor(x))),h,gam,1000)

r <- (log(y) - fit$mu.fitted)^2/fit$xix
plot(x,r,xlim=rx,ylim=range(r),type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(sa[,1],exp(sa[,2]),xlim=rx,ylim=range(r),type="l",ylab="",xlab=xl,main="Skewness")

########################### Residual analysis ##################################

m1 <- "Residuals for the median submodel"
res.dev.mu <- sqrt(fit$deviance.mu)*ifelse(fit$z.hat>=0,1,-1)
yl <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))
plot(x,res.dev.mu,cex=0.3,lwd=3,ylim=yl,main=m1,xlab=xl,ylab="")
abline(h=-3,lty=3)
abline(h=+3,lty=3)
m2 <- "Residuals for the skewness submodel"
res.dev.phi <- sqrt(fit$deviance.phi)*ifelse(fit$z.hat>=0,1,-1)
yl <- c(min(res.dev.phi,-3.5),max(res.dev.phi,3.5))
plot(x,res.dev.phi,cex=0.3,lwd=3,ylim=yl,main=m2,xlab=xl,ylab="")
abline(h=-3,lty=3)
abline(h=+3,lty=3)

########################### Sensitivity analysis ##################################

m1 <- "Local Influence under case-weight perturbation scheme"
m2 <- "Total Local Influence under case-weight perturbation scheme"
plot(fit$cw[,1], type="h", main=m1, xlab="Index", ylab="")
plot(fit$cw[,2], type="h", main=m2, xlab="Index", ylab="")

m1 <- "Local Influence under response perturbation scheme"
m2 <- "Total Local Influence under response perturbation scheme"
plot(fit$pr[,1], type="h", main=m1, xlab="Index", ylab="")
plot(fit$pr[,2], type="h", main=m2, xlab="Index", ylab="")

#########################################################################################
############ European Rabbits Data - a Student-t semiparametric model #############
#########################################################################################

data(Erabbits)
Erabbits2 <- Erabbits[order(Erabbits$age,Erabbits$wlens),]
attach(Erabbits2)

loc.f <- function(vP){
	  beta <- vP[1:3]
	  exp(beta[1] - beta[2]/(beta[3] + age))
}
start <- c(b1=5.6, b2=128, b3=36.45)
fit <- ssym.nl(wlens,loc.f,start=start,ncs=age,family='Student',xi=5,
               local.influence=TRUE)
summary(fit)

####################### Plot of the fitted model ##############################

par(mfrow=c(4,2))
xl <- "Age of the animal (in days)"
rx <- range(age)
ry <- range(wlens)
plot(age,wlens,xlim=rx,ylim=ry,type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(age,fit$mu.fitted,xlim=rx,ylim=ry,type="l",ylab="",xlab=xl,main="Mean")

h <- fit$coefs.phi
ss <- splinek(as.numeric(levels(factor(age))),age)
gam <- solve(ss$R)%*%t(ss$Q)%*%h
sa <- ncs.graph(as.numeric(levels(factor(age))),h,gam,1000)

r <- (wlens - fit$mu.fitted)^2/fit$xix
ry <- range(r)
plot(age,r,xlim=rx,ylim=ry,type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(sa[,1],exp(sa[,2]),xlim=rx,ylim=ry,type="l",ylab="",xlab=xl,main="Dispersion")

########################### Residual analysis ##################################

m1 <- "Residuals for the mean submodel"
res.dev.mu <- sqrt(fit$deviance.mu)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))
plot(age,res.dev.mu,cex=0.3,lwd=3,ylim=ry,main=m1,xlab=xl,ylab="")
abline(h=-3,lty=3)
abline(h=+3,lty=3)
m2 <- "Residuals for the dispersion submodel"
res.dev.phi <- sqrt(fit$deviance.phi)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.phi,-3.5),max(res.dev.phi,3.5))
plot(age,res.dev.phi,cex=0.3,lwd=3,ylim=ry,main=m2,xlab=xl,ylab="")
abline(h=-3,lty=3)
abline(h=+3,lty=3)

########################### Sensitivity analysis ##################################

m1 <- "Local Influence under case-weight perturbation scheme"
m2 <- "Total Local Influence under case-weight perturbation scheme"
plot(fit$cw[,1], type="h", main=m1, xlab="Index", ylab="")
plot(fit$cw[,2], type="h", main=m2, xlab="Index", ylab="")

m1 <- "Local Influence under response perturbation scheme"
m2 <- "Total Local Influence under response perturbation scheme"
plot(fit$pr[,1], type="h", main=m1, xlab="Index", ylab="")
plot(fit$pr[,2], type="h", main=m2, xlab="Index", ylab="")
}