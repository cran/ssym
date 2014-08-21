\name{ssym.l}
\alias{ssym.l}
\title{Fitting Semiparametric Symmetric (or Log-Symmetric) Regression Models}
\description{
\bold{ssym.l} is used to fit a semiparametric regression model suitable for 
data set analysis in which the conditional distribution of the response is 
symmetric and continuous. In this setup, both location and dispersion
parameters of the response variable distribution are explicitly modeled
through semiparametric functions.}
\usage{
ssym.l(response, formula.mu, ncs.mu, start.lambda.mu, lambda.mu, formula.phi, ncs.phi,
start.lambda.phi, lambda.phi, family, xi, epsilon, maxiter, subset, local.influence)
}
\arguments{
  \item{response}{the response (or log-response) variable.}
  \item{formula.mu}{(optional) a symbolic description of the parametric function to be fitted to the location parameter.}
  \item{ncs.mu}{(optional) an explanatory (continuous) variable to be used in the nonparametric function to be fitted to the location parameter.}
  \item{start.lambda.mu}{(optional) a numeric value of starting estimate for the smoothing parameter associated to the location submodel.}
  \item{lambda.mu}{(optional) a numerical value for the smoothing parameter indicating that it is provided by the user rather than estimated from the data.}
  \item{formula.phi}{(optional) a symbolic description of the parametric function to be fitted to the dispersion parameter.}
  \item{ncs.phi}{(optional) an explanatory (continuous) variable to be used in the nonparametric function to be fitted to the dispersion parameter.}
  \item{start.lambda.phi}{(optional) a numeric value of starting estimate for the smoothing parameter associated to the location submodel.}
  \item{lambda.phi}{(optional) a numerical value for the smoothing parameter indicating that it is provided by the user rather than estimated from the data.}
  \item{family}{a description of the error distribution to be used in the model. Supported families include \emph{Normal}, \emph{Student}, \emph{Contnormal},
               \emph{Powerexp}, \emph{Hyperbolic}, \emph{Slash}, \emph{Sinh-normal} and \emph{Sinh-t}, which correspond to normal, Student-t, contaminated normal,
			    power exponential, symmetric hyperbolic, slash, sinh-normal and sinh-t distributions, respectively.}
  \item{xi}{a numeric value or numeric vector that represents the extra parameter of the specified error distribution.}
  \item{epsilon}{(optional) positive convergence tolerance. Default value is 0.0000001.}
  \item{maxiter}{(optional) an positive integer giving the maximal number of iterations for the estimating process. Default value is 1000.}
  \item{subset}{(optional) expression indicating individuals or observations to keep (or drop).}
  \item{local.influence}{(optional) logical. If TRUE, local influence measures for the location parameters vector under two perturbation schemes are calculated.}
}
\details{
The iterative estimation process for the interest parameters is based on the Fisher scoring and backfitting algorithms. Because some distributions such as
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
\item{lambda.mu}{estimate for the smoothing parameter associated to the location parameter.}
\item{dfe.mu}{degrees of freedom associated with the nonparametric part of the location submodel.}
\item{lambda.phi}{estimate for the smoothing parameter associated to the dispersion parameter.}
\item{dfe.phi}{degrees of freedom associated with the nonparametric part of the dispersion submodel.}
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
\code{\link{ssym.nl}}
}

\examples{
###################################################################################
########### Fraction of Cell Volume Data - a Power-exponential model  ###########
###################################################################################

data(Ovocytes)
attach(Ovocytes)

fit <- ssym.l(fraction, formula.mu=~type, ncs.mu=time, formula.phi=~type,
              family="Powerexp", xi=-0.7, local.influence=TRUE)
summary(fit)

####################### Plot of the fitted model ##############################

par(mfrow=c(1,2))
xl <- "Time"
rx <- range(time)
ry <- range(fraction)
h <- fit$coefs.mu[2:length(fit$coefs.mu)]
sa <- ncs.graph(time,h,1000)

plot(time[type=='Mature'],fraction[type=='Mature'],xlim=rx,ylim=ry,type="p",
     cex=0.5,lwd=1,ylab="",xlab="")
par(new=TRUE)
plot(time[type=='Immature'],fraction[type=='Immature'],xlim=rx,ylim=ry,
     type="p",cex=0.5,lwd=2,ylab="",xlab="")
par(new=TRUE)
plot(sa[,1],(sa[,2]),xlim=rx,ylim=ry,type="l",ylab="",xlab=xl,lty=1,main="")
par(new=TRUE)
plot(sa[,1],(sa[,2] + fit$coefs.mu[1]),xlim=rx,ylim=ry,type="l",ylab="",xlab=xl,
     lty=1,main="Mean")
legend(rx[1],ry[2],pt.lwd=c(1,2),bty="n",legend=c("Mature","Immature"),
       pt.cex=0.5,pch=1)

r <- (fraction - fit$mu.fitted)^2/fit$xix
ry <- range(r)
plot(time[type=='Mature'],r[type=='Mature'],xlim=rx,ylim=ry,type="p",cex=0.5,
     lwd=1,ylab="",xlab="")
par(new=TRUE)
plot(time[type=='Immature'],r[type=='Immature'],xlim=rx,ylim=ry,type="p",
     cex=0.5,lwd=2,ylab="",xlab="",main="Dispersion")
abline(h=exp(fit$coefs.phi[1]))
abline(h=exp(fit$coefs.phi[1] + fit$coefs.phi[2]))
legend(rx[1],ry[2],pt.lwd=c(1,2),bty="n",legend=c("Mature","Immature"),
       pt.cex=0.5,pch=1)

################# Weights and Overall goodness-of-fit #########################

par(mfrow=c(1,2))
plot(fit$z.hat,fit$weights,cex=0.3,lwd=3,main="Weights")
d1 <- qqnorm(qnorm(fit$cdfz),plot.it=TRUE,main="Overall Goodness-of-fit measure")
abline(0,1,lty=3)
	   
########################### Residual analysis ##################################

par(mfrow=c(1,2))
m1 <- "Residuals for the mean submodel"
res.dev.mu <- sqrt(fit$deviance.mu)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))
plot(time,res.dev.mu,cex=0.3,lwd=3,ylim=ry,main=m1,xlab=xl,ylab="")
abline(h=-3,lty=3)
abline(h=+3,lty=3)
m2 <- "Residuals for the dispersion submodel"
res.dev.phi <- sqrt(fit$deviance.phi)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.phi,-3.5),max(res.dev.phi,3.5))
plot(time,res.dev.phi,cex=0.3,lwd=3,ylim=ry,main=m2,xlab=xl,ylab="")
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

#########################################################################################
################## Biaxial Fatigue Data - a Birnbaum-Saunders model   ##################
#########################################################################################

data(Biaxial)
attach(Biaxial)

fit <- ssym.l(log(Life), formula.mu=~log(Work), family='Sinh-normal', xi=1.54,
              local.influence=TRUE)
summary(fit)

####################### Plot of the fitted model ##############################

xl <- "Work per cycle"
rx <-range(Work)
ry <- range(Life)
plot(Work,Life,xlim=rx,ylim=ry,type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(Work,exp(fit$mu.fitted),xlim=rx,ylim=ry,type="l",ylab="Life",xlab=xl,main="Median")

################# Weights and Overall goodness-of-fit #########################

par(mfrow=c(1,2))
plot(fit$z.hat,fit$weights,cex=0.3,lwd=3,main="Weights")
d1 <- qqnorm(qnorm(fit$cdfz),plot.it=TRUE,main="Overall Goodness-of-fit measure")
abline(0,1,lty=3)

########################### Residual analysis ##################################

res.dev.mu <- sqrt(fit$deviance.mu)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))
plot(Work,res.dev.mu,cex=0.3,lwd=3,ylim=ry,main="Residuals",xlab=xl,ylab="")
abline(h=-3,lty=3)
abline(h=+3,lty=3)

########################### Sensitivity analysis ##################################

par(mfrow=c(1,2))
m1 <- "Local Influence under case-weight perturbation scheme"
m2 <- "Total Local Influence under case-weight perturbation scheme"
plot(fit$cw[,1], type="h", main=m1, xlab="Index", ylab="")
plot(fit$cw[,2], type="h", main=m2, xlab="Index", ylab="")
}
