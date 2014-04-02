\name{ssym.l}
\alias{ssym.l}
\title{Fitting Semiparametric Symmetric (or Log-Symmetric) Regression Models}
\description{
\bold{ssym.l} is used to fit a semiparametric regression model suitable for 
data set analysis in which the conditional distribution of the response (or log-response) is 
symmetric and continuous. In this setup, both location and dispersion
parameters of response (or log-response) variable distribution are explicitly modeled, 
the location using a linear function and the dispersion using a 
semiparametric function, which is aproximated by a natural cubic spline.}
\usage{
ssym.l(response, formula.mu, start, formula.phi, ncs, start.lambda, lambda, family, xi,
       epsilon, maxiter, subset, local.influence)
}
\arguments{
  \item{response}{the response (or log-response) variable.}
  \item{formula.mu}{a symbolic description of the parametric function to be fitted to the location parameter.}
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
\code{\link{ssym.nl}}
}

\examples{
#########################################################################################
############## Texture of Snacks Data - a log-normal semiparametric model ###############
#########################################################################################

data(Snacks)
Snacks2 <- Snacks[order(Snacks$snack,Snacks$week,Snacks$texture),]
attach(Snacks2)

snack <- factor(snack)
week2 <- week^2
week3 <- week^3

fit <- ssym.l(log(texture),formula.mu=~snack + week2 + week3,formula.phi=~snack,
              ncs=week,family='Normal',local.influence=TRUE)
summary(fit)

####################### Plot of the fitted model ##############################

par(mfrow=c(4,2))
rx <- range(week)
ry <- range(texture)
plot(week,texture,xlim=rx,ylim=ry,type="p",cex=0.3,lwd=3,ylab="",xlab="")
for(i in 1:5){
 par(new=TRUE)
 plot(week[snack==i],exp(fit$mu.fitted[snack==i]),xlim=rx,ylim=ry,
      type="l",ylab="",xlab="",lty=i,main="Median")
}
legend(min(week),max(texture),lty=1:5,bty="n",legend=paste("snack",1:5))

h <- fit$coefs.phi[5:length(fit$coefs.phi)]
ss <- splinek(as.numeric(levels(factor(week))),week)
gam <- solve(ss$R)%*%t(ss$Q)%*%h
sa <- ncs.graph(as.numeric(levels(factor(week))),h,gam,1000)

r <- (log(texture) - fit$mu.fitted)^2/fit$xix
ry <- range(r)
plot(week,r,xlim=rx,ylim=ry,type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(sa[,1],exp(sa[,2]),xlim=rx,ylim=ry,type="l",ylab="",xlab="",lty=1)
for(i in 1:4){
 par(new=TRUE)
 plot(sa[,1],exp(sa[,2]+fit$coefs.phi[i]),xlim=rx,ylim=ry,type="l",ylab="",xlab="",
      lty=i,main="Skewness")
}
legend(min(week),max(r),lty=1:5,bty="n",legend=paste("snack",1:5))

########################### Residual analysis ##################################

xl <- "Week"
m1 <- "Residuals for the median submodel"
res.dev.mu <- sqrt(fit$deviance.mu)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))
plot(week,res.dev.mu,cex=0.3,lwd=3,ylim=ry,main=m1,xlab=xl,ylab="")
abline(h=-3,lty=3)
abline(h=+3,lty=3)
m2 <- "Residuals for the skewness submodel"
res.dev.phi <- sqrt(fit$deviance.phi)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.phi,-3.5),max(res.dev.phi,3.5))
plot(week,res.dev.phi,cex=0.3,lwd=3,ylim=ry,main=m2,xlab=xl,ylab="")
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
################## Biaxial Fatigue Data - a Birnbaum-Saunders model   ##################
#########################################################################################

data(Biaxial)
Biaxial2 <- Biaxial[order(Biaxial$Work,Biaxial$Life),]
attach(Biaxial2)

fit <- ssym.l(log(Life), formula.mu=~log(Work), family='Sinh-normal', xi=1.54,
              local.influence=TRUE)
summary(fit)

####################### Plot of the fitted model ##############################

par(mfrow=c(3,2))
xl <- "Work per cycle"
rx <-range(Work)
ry <- range(Life)
plot(Work,Life,xlim=rx,ylim=ry,type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
plot(Work,exp(fit$mu.fitted),xlim=rx,ylim=ry,type="l",ylab="Life",xlab=xl,main="Median")

########################### Residual analysis ##################################

res.dev.mu <- sqrt(fit$deviance.mu)*ifelse(fit$z.hat>=0,1,-1)
ry <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))
plot(Work,res.dev.mu,cex=0.3,lwd=3,ylim=ry,main="Residuals",xlab=xl,ylab="")
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
