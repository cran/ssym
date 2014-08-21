\name{ssym.nl}
\alias{ssym.nl}
\title{Fitting Semiparametric Symmetric (or Log-Symmetric) Regression Models}
\description{
\bold{ssym.nl} is used to fit a semiparametric regression model suitable for 
data set analysis in which the conditional distribution of the response is 
symmetric and continuous. In this setup, both location and dispersion
parameters of the response variable distribution are explicitly modeled, 
the location using a nonlinear function and the dispersion using a 
semiparametric function.}
\usage{
ssym.nl(response, mu, start, formula.phi, ncs.phi, start.lambda.phi, lambda.phi,
        family, xi,	epsilon, maxiter, subset, local.influence)
}
\arguments{
  \item{response}{the response (or log-response) variable.}
  \item{mu}{an R function with the nonlinear model that describes the location parameter of the response (or log-response) variable distribution.}
  \item{start}{a named numeric vector of starting estimates.}
  \item{formula.phi}{(optional) a symbolic description of the parametric function to be fitted to the dispersion parameter.}
  \item{ncs.phi}{(optional) an explanatory (continuous) variable to be used in the nonparametric function to be fitted to the dispersion parameter.}
  \item{start.lambda.phi}{(optional) a numeric value of starting estimate for the smoothing parameter.}
  \item{lambda.phi}{(optional) a numerical value for the smoothing parameter indicating that it is provided by the user rather than estimated from the data.}
  \item{family}{a description of the error distribution to be used in the model. Supported families include \emph{Normal}, \emph{Student}, \emph{Contnormal},
               \emph{Powerexp}, \emph{Hyperbolic}, \emph{Slash}, \emph{Sinh-normal} and \emph{Sinh-t}, which correspond to normal, Student-t, contaminated normal,
			    power exponential, symmetric hyperbolic, slash, sinh-normal and sinh-t distributions, respectively.}
  \item{xi}{a numeric value or numeric vector that represents the extra parameter of the specified error distribution.}
  \item{epsilon}{(optional) positive convergence tolerance. Default value is 0.0000001.}
  \item{maxiter}{(optional) an positive integer giving the maximal number of iterations for the estimating process. Default value is 500.}
  \item{subset}{(optional) expression indicating individuals or observations to keep (or drop).}
  \item{local.influence}{(optional) logical. If TRUE, local influence measures for the location parameters vector under two perturbation schemes are calculated.}
}
\details{
The iterative estimation process for the parameters of interest is based on the Fisher scoring and backfitting algorithms. Because some distributions such as
Student-t, contaminated normal, power exponential, slash and symmetric hyperbolic may be obtained as a scale mixture of normal distributions, the EM algorithm
is applied in those cases to obtain a more efficient iterative process for the parameter estimation. Further, because the Sinh-t distribution can be obtained as a 
scale mixture of Sinh-normal distributions, the EM algorithm is also applied in that case to obtain a more efficient iterative process for the parameter estimation.
The smoothing parameter is chosen using the unweighted cross-validation score.
}
\value{
\item{coefs.mu}{a vector of parameter estimates associated with the nonlinear function fitted to the location of the response (or log-response) variable.}
\item{coefs.phi}{a vector of parameter estimates associated with the semiparametric function fitted to the dispersion of the response (or log-response) variable.}
\item{se.mu}{approximated standard errors asociated with \emph{coefs.mu}.}
\item{se.phi}{approximated standard errors asociated with \emph{coefs.phi}.}
\item{weights}{final weights of the iterative process.}
\item{lambda.phi}{estimate for the smoothing parameter.}
\item{dfe.phi}{degrees of freedom associated with the nonparametric part of the model.}
\item{deviance.mu}{a vector of \emph{deviances} associated with the location of the response (or log-response) variable.}
\item{deviance.phi}{a vector of \emph{deviances} associated with the dispersion of the response (or log-response) variable.}
\item{mu.fitted}{a vector of fitted values for the location of the response (or log-response) variable.}
\item{phi.fitted}{a vector of fitted values for the dispersion of the response (or log-response) variable.}
\item{lpdf}{a vector of individual contributions to the log-likelihood function.}
\item{cdfz}{a vector of the cumulative distribution function of each individual.}
\item{cw}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the case-weight perturbation scheme)
                   associated with \emph{coefs.mu} and based on the conformal normal curvature.}
\item{pr}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the response perturbation scheme)
                   associated with \emph{coefs.mu} and based on the conformal normal curvature.}
\item{cw.theta}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the case-weight perturbation scheme)
                asociated with the (entire) vector of parameters and based on the conformal normal curvature.}
\item{pr.theta}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the response perturbation scheme)
                asociated with the (entire) vector of parameters and based on the conformal normal curvature.}
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

par(mfrow=c(1,2))
xl <- "Distance to the metal"
rx <- range(x)
plot(x,y,xlim=rx,ylim=range(y),type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(x,exp(fit$mu.fitted),xlim=rx,ylim=range(y),type="l",ylab="",xlab=xl,main="Median")

h <- fit$coefs.phi
sa <- ncs.graph(x,h,1000)

r <- (log(y) - fit$mu.fitted)^2/fit$xix
plot(x,r,xlim=rx,ylim=range(r),type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(sa[,1],exp(sa[,2]),xlim=rx,ylim=range(r),type="l",ylab="",xlab=xl,main="Skewness")

################# Weights and Overall goodness-of-fit #########################

par(mfrow=c(1,2))
plot(fit$z.hat,fit$weights,cex=0.3,lwd=3,main="Weights")
d1 <- qqnorm(qnorm(fit$cdfz),plot.it=TRUE,main="Overall Goodness-of-fit measure")
abline(0,1,lty=3)

########################### Residual analysis ##################################

par(mfrow=c(1,2))
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

par(mfrow=c(1,2))
m1 <- "Local Influence under case-weight perturbation scheme"
m2 <- "Total Local Influence under case-weight perturbation scheme"
plot(fit$cw[,1], type="h", main=m1, xlab="Index", ylab="")
plot(fit$cw[,2], type="h", main=m2, xlab="Index", ylab="")

par(mfrow=c(1,2))
m1 <- "Local Influence under response perturbation scheme"
m2 <- "Total Local Influence under response perturbation scheme"
plot(fit$pr[,1], type="h", main=m1, xlab="Index", ylab="")
plot(fit$pr[,2], type="h", main=m2, xlab="Index", ylab="")

########################################################################################
############### Body Fat percentage Data - a Birnbaum-Saunders-t model #################
#########################################################################################

library(sn)
data(ais)
y <- sort(ais$Bfat)

location.f <- function(vP){rep(vP[1],length(y))}

fit <- ssym.nl(log(y), mu=location.f, start=mean(log(y)), family="Sinh-t", xi=c(4.5,4))
summary(fit)

xl <- "Body Fat percentage"
ss <- c(seq(5,10,by=2.5),seq(15,20,by=5),seq(22.5,37.5,by=2.5))
par(mfrow=c(1,2))
hist(y,xlim=range(y),ylim=c(0,0.1),prob=TRUE,breaks=ss,col="light gray",
     border="dark gray",xlab="",main="",ylab="")
par(new=TRUE)
plot(y,exp(fit$lpdf)/y,xlim=range(y),ylim=c(0,0.1),type="l",
     xlab=xl,ylab="",main="Histogram")

plot(y,fit$cdfz,xlim=range(y),ylim=c(0,1),type="l",xlab="",ylab="")
par(new=TRUE)
plot(ecdf(y),xlim=range(y),ylim=c(0,1),verticals=TRUE,do.points=FALSE,
     col="dark gray",xlab=xl,main="ecdf")

########################################################################################
######### Gross Domestic Product per capita Data - a Birnbaum-Saunders model ###########
#########################################################################################

data(gdp)
attach(gdp)

location.f <- function(vP){rep(vP[1],length(gdp2010))}

fit <- ssym.nl(log(gdp2010), mu=location.f, start=mean(log(gdp2010)), family='Sinh-normal', xi=2.2)
summary(fit)

xl <- "GDP per capita 2010"
par(mfrow=c(1,2))
hist(gdp2010,xlim=range(gdp2010),ylim=c(0,0.00015),prob=TRUE,breaks=55,col="light gray",
     border="dark gray",xlab="",main="",ylab="")
par(new=TRUE)
plot(gdp2010,exp(fit$lpdf)/gdp2010,xlim=range(gdp2010),ylim=c(0,0.00015),type="l",
     xlab=xl,ylab="",main="Histogram")

plot(gdp2010,fit$cdfz,xlim=range(gdp2010),ylim=c(0,1),type="l",xlab="",ylab="")
par(new=TRUE)
plot(ecdf(gdp2010),xlim=range(gdp2010),ylim=c(0,1),verticals=TRUE,do.points=FALSE,
     col="dark gray",xlab=xl,main="ecdf")
}