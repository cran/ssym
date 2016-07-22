\name{ssym.l}
\alias{ssym.l}
\title{Fitting Semi-parametric Log-symmetric Regression Models}
\description{
\bold{ssym.l} is used to fit a semi-parametric regression model suitable for 
analysis of data sets in which the response variable is continuous, strictly
positive, and asymmetric. Under this setup, both median and skewness of the response
variable distribution are explicitly modeled through semi-parametric functions,
whose nonparametric components may be approximated by natural cubic splines or P-splines.}
\usage{
ssym.l(formula, family, xi, data, epsilon, maxiter, subset, link.mu, link.phi,
       local.influence, spec, std.out)
}
\arguments{
  \item{formula}{a symbolic description of the systematic component of the model to be fitted. See details for further information.}
  \item{family}{a description of the (log) error distribution to be used in the model. Supported families include \code{Normal}, \code{Student}, \code{Contnormal},
               \code{Powerexp}, \code{Hyperbolic}, \code{Slash}, \code{Sinh-normal} and \code{Sinh-t}, which correspond to normal, Student-t, contaminated normal,
			    power exponential, symmetric hyperbolic, slash, sinh-normal and sinh-t distributions, respectively.}
  \item{xi}{a numeric value or numeric vector that represents the extra parameter of the specified error distribution.}
  \item{data}{an optional data frame, list or environment containing the variables in the model.}
  \item{epsilon}{an optional positive value, which represents the convergence criterion. Default value is 1e-07.}
  \item{maxiter}{an optional positive integer giving the maximal number of iterations for the estimating process. Default value is 1e03.}
  \item{subset}{an optional expression that specifies a subset of individuals to be used in the fitting process.}
  \item{link.mu}{an optional character that specifies the link function of the median submodel.}
  \item{link.phi}{an optional character that specifies the link function of the skewness submodel.}
  \item{local.influence}{logical. If TRUE, local influence measures under two perturbation schemes are calculated. Default is \code{FALSE}.}
  \item{spec}{an optional character. The smoothing parameter is estimated by minimizing a overall goodness-of-fit criterion such as AIC or BIC. \code{spec} is an optional
              string to specify the goodness-of-fit measure to be used. Default value is \code{AIC}.}
  \item{std.out}{logical. If FALSE, just a reduced set of attributes is returned by the model-fitting function. Default is \code{TRUE}}.
}
\details{
The argument \emph{formula} comprises of three parts (separated by the symbols "~" and "|"), namely: observed response variable in log-scale, predictor of the
median submodel (having logarithmic link) and predictor of the skewness (or the relative dispersion) submodel (having logarithmic link). An arbitrary number of
nonparametric effects may be specified in the predictors. These effects are specified to be approximated by natural cubic splines or P-splines using the functions
\code{ncs()} or \code{psp()}, respectively.

The iterative estimation process is based on the Fisher scoring and backfitting algorithms. Because some distributions such as log-Student-t,
log-contaminated-normal, log-power-exponential, log-slash and log-hyperbolic may be obtained as a power mixture of the log-normal distribution, the
expectation-maximization (EM) algorithm is applied in those cases to obtain a more efficient iterative process of parameter estimation. Furthermore,
because the Birnbaum-Saunders-t distribution can be obtained as a scale mixture of the Birnbaum-Saunders distribution, the expectation-maximization algorithm
is also applied in this case to obtain a more efficient iterative process of parameter estimation. The smoothing parameter is chosen by minimizing
the AIC or BIC criteria.

The function \code{ssym.l()} calculates overall goodness-of-fit statistics, deviance-type residuals for both submodels, as well as local influence measures
under the case-weight and response perturbation schemes.
}

\value{
\item{theta.mu}{a vector of parameter estimates associated with the median submodel.}
\item{theta.phi}{a vector of parameter estimates associated with the skewness (or the relative dispersion) submodel.}
\item{vcov.mu}{approximate variance-covariance matrix associated with the median submodel.}
\item{vcov.phi}{approximate variance-covariance matrix associated with the skewness (or the relative dispersion) submodel.}
\item{weights}{final weights of the iterative process.}
\item{lambdas.mu}{estimate of the smoothing parameter(s) associated with the nonparametric part of the median submodel.}
\item{lambdas.phi}{estimate of the smoothing parameter(s) associated with the nonparametric part of the skewness (or the relative dispersion) submodel.}
\item{gle.mu}{degrees of freedom associated with the nonparametric part of the median submodel.}
\item{gle.phi}{degrees of freedom associated with the nonparametric part of the skewness (or the relative dispersion) submodel.}
\item{deviance.mu}{a vector with the individual contributions to the deviance associated with the median submodel.}
\item{deviance.phi}{a vector with the individual contributions to the deviance associated with the skewness (or the relative dispersion) submodel.}
\item{mu.fitted}{a vector with the fitted values of the (in log-scale) median submodel.}
\item{phi.fitted}{a vector with the fitted values of the skewness (or the relative dispersion) submodel.}
\item{lpdf}{a vector of individual contributions to the log-likelihood function.}
}

\references{Vanegas, L.H. and Paula, G.A. (2015) A semiparametric approach for joint modeling of median and skewness. TEST 24, 110-135.

			Vanegas, L.H. and Paula, G.A. (2016) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics 30, 196-220.

			Vanegas, L.H. and Paula, G.A. (2016) An extension of log-symmetric regression models: R codes and applications. Journal of Statistical Computation and Simulation 86, 1709-1735.}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\seealso{
\code{\link{ssym.nl}, \link{ssym.l2}}
}

\examples{

###################################################################################
######### Fraction of Cell Volume Data - a log-power-exponential model  ###########
###################################################################################
#data("Ovocytes", package="ssym")
#fit <- ssym.l(log(fraction) ~ type + psp(time) | type + psp(time),
#               data=Ovocytes, family='Powerexp', xi=-0.55, local.influence=TRUE)
#summary(fit)
#
################## Graph of the nonparametric effects ##################
#par(mfrow=c(1,2))
#np.graph(fit, which=1, exp=TRUE)
#np.graph(fit, which=2, exp=TRUE)
#
################## Graph of deviance-type residuals ##################
#plot(fit)
#
################### Simulated envelopes ##################
#envelope(fit)
#
################### Graph of local influence measures ##################
#out <- influence(fit)

###################################################################################
############### Textures of snacks Data - a log-Student-t model  #################
###################################################################################
#data("Snacks", package="ssym")
#fit <- ssym.l(log(texture) ~ type + ncs(week) | type, data=Snacks,
#        family='Student', xi=15, local.influence=TRUE)
#summary(fit)
#
############################ Extra parameter ###########################
#extra.parameter(fit,5,50)
#
################### Graph of the nonparametric effect ##################
#np.graph(fit, which=1, exp=TRUE)
#
################### Graph of deviance-type residuals ##################
#plot(fit)
#
################### Simulated envelopes ##################
#envelope(fit)
#
################### Plot of influence measures ##################
#out <- influence(fit)

###################################################################################
################### Daphnia Data - a log-normal model ########################
###################################################################################
#data("daphnia", package="nlreg")
#fit <- ssym.l(log(time) ~ ncs(conc) | ncs(conc), data=daphnia, family="Normal")
#summary(fit)
#
################### Graph of the nonparametric effects ##################
#par(mfrow=c(1,2))
#np.graph(fit, which=1, exp=TRUE)
#np.graph(fit, which=2, exp=TRUE)
#
################### Simulated envelopes ##################
#envelope(fit)

###################################################################################
####################### gam.data - a Power-exponential model   ####################
###################################################################################
#data("gam.data", package="gam")
#
#fit <- ssym.l(y~psp(x),data=gam.data,family="Powerexp",xi=-0.5)
#summary(fit)
#
################## Graph of the nonparametric effect ##################
#np.graph(fit, which=1)
#
###################################################################################
######### Personal Injury Insurance Data - a Birnbaum-Saunders-t model   ##########
###################################################################################
#data("Claims", package="ssym")
#fit <- ssym.l(log(total) ~ op_time | op_time, data=Claims,
#        family='Sinh-t', xi=c(0.1,4), local.influence=TRUE)
#summary(fit)
#
################## Plot of deviance-type residuals ##################
#plot(fit)
#
################### Simulated envelopes ##################
#envelope(fit)
################## Plot of influence measures ##################
#out <- influence(fit)

###################################################################################
######### Body Fat Percentage Data - a Birnbaum-Saunders-t model   ##########
###################################################################################
#data("ais", package="sn")
#fit <- ssym.l(log(Bfat)~1, data=ais, family='Sinh-t', xi=c(4.5,4))
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,c(3,4),c(5,7))
#
################## Plot of the fitted model ##################
#id <- sort(ais$Bfat, index=TRUE)$ix
#par(mfrow=c(1,2))
#hist(ais$Bfat[id],xlim=range(ais$Bfat),ylim=c(0,0.1),prob=TRUE,breaks=15,
#     col="light gray",border="dark gray",xlab="",ylab="",main="")
#par(new=TRUE)
#plot(ais$Bfat[id],exp(fit$lpdf[id])/ais$Bfat[id],xlim=range(ais$Bfat),
#     ylim=c(0,0.1),type="l",xlab="",ylab="Density",main="Histogram")
#	 
#plot(ais$Bfat[id],fit$cdfz[id],xlim=range(ais$Bfat),ylim=c(0,1),type="l",
#     xlab="",ylab="",main="")
#par(new=TRUE)
#plot(ecdf(ais$Bfat[id]),xlim=range(ais$Bfat),ylim=c(0,1),verticals=TRUE,
#     do.points=FALSE,col="dark gray",ylab="Probability",xlab="",main="ECDF")

###################################################################################
################### ALCOA Aluminium Data - a log-slash model   ####################
###################################################################################

#data("alcoa", package="robustloggamma")
#alcoa2 <- data.frame(alcoa$dist[alcoa$label=="C"])
#colnames(alcoa2) <- "dist"
#
#fit <- ssym.l(log(dist) ~ 1, data=alcoa2, family="Slash", xi=1.212)
#
################## Plot of the fitted model ##################
#id <- sort(alcoa2$dist, index=TRUE)$ix
#par(mfrow=c(1,2))
#hist(alcoa2$dist[id],xlim=c(0,45),ylim=c(0,0.1),prob=TRUE,breaks=60,
#     col="light gray",border="dark gray",xlab="",ylab="",main="")
#par(new=TRUE)
#plot(alcoa2$dist[id],exp(fit$lpdf[id])/alcoa2$dist[id],xlim=c(0,45),
#ylim=c(0,0.1), type="l",xlab="",ylab="",main="")
#	 
#plot(alcoa2$dist[id],fit$cdfz[id],xlim=range(alcoa2$dist),ylim=c(0,1),type="l",
#     xlab="",ylab="",main="")
#par(new=TRUE)
#plot(ecdf(alcoa2$dist[id]),xlim=range(alcoa2$dist),ylim=c(0,1),verticals=TRUE,
#     do.points=FALSE,col="dark gray",ylab="",xlab="",main="")

##################################################################################
################### Boston Housing Data - a log-Slash model   ####################
###################################################################################
#data("Boston", package="MASS")
#fit <- ssym.l(log(medv) ~ crim + rm + tax + psp(lstat) + psp(dis) | psp(lstat),
#              data=Boston, family="Slash", xi=1.56, local.influence=TRUE)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,1.0,2.3)
#
################## Plot of deviance-type residuals ##################
#plot(fit)
#
################## Plot of nonparametric effects ##################
#par(mfrow=c(1,3))
#np.graph(fit,which=1,exp=TRUE,"lstat")
#np.graph(fit,which=1,exp=TRUE,"dis")
#np.graph(fit,which=2,exp=TRUE,"lstat")
#
################## Plot of influence measures ##################
#out <- influence(fit)
#
################### Simulated envelopes ##################
#envelope(fit)

###################################################################################
####################### mcycle Data - a Power-exponential model   #################
###################################################################################
#data("mcycle", package="MASS")
#fit <- ssym.l(accel ~ ncs(times)|ncs(times), data=mcycle, family="Powerexp",xi=-0.6)
#summary(fit)
#
################## Plot of nonparametric effects ##################
#par(mfrow=c(1,2))
#np.graph(fit,which=1,obs=TRUE)
#np.graph(fit,which=2,exp=TRUE,obs=TRUE)
#
################### Simulated envelopes ##################
#envelope(fit)

###################################################################################
################### Steel Data - a log-hyperbolic model   ####################
###################################################################################
#data("Steel", package="ssym")
#fit <- ssym.l(log(life)~psp(stress), data=Steel, family="Hyperbolic", xi=1.25)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,0.5,2)
#
################## Plot of nonparametric effects ##################
#np.graph(fit,which=1,exp=TRUE)
	 
}

