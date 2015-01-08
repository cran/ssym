\name{ssym.l}
\alias{ssym.l}
\title{Fitting Semi-parametric Log-symmetric Regression Models}
\description{
\bold{ssym.l} is used to fit a semi-parametric regression model suitable for 
analysis of data sets in which the response variable is continuous, strictly
positive, and asymmetric. In this setup, both median and skewness of the response
variable distribution are explicitly modeled through semi-parametric functions,
whose nonparametric components may be approximated by natural cubic splines or P-splines.}
\usage{
ssym.l(formula, family, xi, data, epsilon, maxiter, subset, local.influence)
}
\arguments{
  \item{formula}{a symbolic description of the systematic component of the model to be fitted. See details for further information.}
  \item{family}{a description of the (log) error distribution to be used in the model. Supported families include \emph{Normal}, \emph{Student}, \emph{Contnormal},
               \emph{Powerexp}, \emph{Hyperbolic}, \emph{Slash}, \emph{Sinh-normal} and \emph{Sinh-t}, which correspond to normal, Student-t, contaminated normal,
			    power exponential, symmetric hyperbolic, slash, sinh-normal and sinh-t distributions, respectively.}
  \item{xi}{a numeric value or numeric vector that represents the extra parameter of the specified error distribution.}
  \item{data}{an optional data frame, list or environment containing the variables in the model.}
  \item{epsilon}{an optional positive value, which represents the convergence criterion. Default value is 1e-07.}
  \item{maxiter}{an optional positive integer giving the maximal number of iterations for the estimating process. Default value is 1e03.}
  \item{subset}{an optional expression specifying a subset of individuals to be used in the fitting process.}
  \item{local.influence}{logical. If TRUE, local influence measures under two perturbation schemes are calculated.}
}
\details{
The argument \emph{formula} comprises of three parts (separated by the symbols "~" and "|"), namely: observed response variable in log-scale, predictor of the
median submodel (having logarithmic link) and predictor of the skewness submodel (having logarithmic link). A non-parametric effect may be
specified in the predictors, either approximated by a natural cubic spline or a P-spline using the functions \emph{ncs()} or \emph{psp()},
respectively.

The iterative estimation process is based on the Fisher scoring and backfitting algorithms. Because some distributions such as log-Student-t,
log-contaminated-normal, log-power-exponential, log-slash and log-hyperbolic may be obtained as a power mixture of the log-normal distribution, the
expectation-maximization (EM) algorithm is applied in those cases to obtain a more efficient iterative process for the parameter estimation. Furthermore,
because the Birnbaum-Saunders-t distribution can be obtained as a scale mixture of the Birnbaum-Saunders distribution, the expectation-maximization algorithm
is also applied in this case to obtain a more efficient iterative process for the parameter estimation. The smoothing parameter(s) is(are)
chosen using the unweighted cross-validation score.

The function \emph{ssym.l()} calculates deviance-type residuals for both submodels as well as local influence measures under case-weight and
response perturbation schemes.
}

\value{
\item{coefs.mu}{a vector of parameter estimates associated with the median submodel.}
\item{coefs.phi}{a vector of parameter estimates associated with the skewness submodel.}
\item{vcov.mu}{approximate variance-covariance matrix associated with the median submodel.}
\item{vcov.phi}{approximate variance-covariance matrix associated with the skewness submodel.}
\item{weights}{final weights of the iterative process.}
\item{lambda.mu}{estimate for the smoothing parameter associated with the nonparametric part of the median submodel.}
\item{dfe.mu}{degrees of freedom associated with the nonparametric part of the median submodel.}
\item{lambda.phi}{estimate for the smoothing parameter associated  with the nonparametric part of the skewness submodel.}
\item{dfe.phi}{degrees of freedom associated with the nonparametric part of the skewness submodel.}
\item{deviance.mu}{a vector of \emph{deviances} associated with the median submodel.}
\item{deviance.phi}{a vector of \emph{deviances} associated with the skewness submodel.}
\item{mu.fitted}{a vector of fitted values for the (log) median submodel.}
\item{phi.fitted}{a vector of fitted values for the skewness submodel.}
\item{lpdf}{a vector of individual contributions to the log-likelihood function.}

\item{cw}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the case-weight perturbation scheme)
                   associated with the median submodel.}
\item{pr}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the response perturbation scheme)
                   associated with the median submodel.}
\item{cw.theta}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the case-weight perturbation scheme).}
\item{pr.theta}{if \emph{local.influence=TRUE}, a matrix of local influence and total local influence measures (under the response perturbation scheme).}
}

\references{Vanegas, L.H. and Paula, G.A. (2015a) A Semiparametric Approach for Joint Modeling of Median and Skewness. TEST (to appear)


			Vanegas, L.H. and Paula, G.A. (2015b) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics (to appear)}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\seealso{
\code{\link{ssym.nl}}
}

\examples{

###################################################################################
######### Fraction of Cell Volume Data - a log-power-exponential model  ###########
###################################################################################

data("Ovocytes", package="ssym")
fit <- ssym.l(log(fraction) ~ type + psp(time) | type + psp(time), data=Ovocytes,
              family='Powerexp', xi=-0.55, maxiter=5000, local.influence=TRUE)
summary(fit)

################## Graph of the nonparametric effects ##################

par(mfrow=c(1,2))
np.graph(fit, which=1, exp=TRUE)
np.graph(fit, which=2, exp=TRUE)

################## Graph of deviance-type residuals ##################

plot(fit)

################## Graph of local influence measures ##################

ilm <- influence.ssym(fit)

###################################################################################
############### Textures of snacks Data - a log-Student-t model  #################
###################################################################################

data("Snacks", package="ssym")
fit <- ssym.l(log(texture) ~ type + ncs(week) | type, data=Snacks,
              family='Student', xi=15)
summary(fit)

################## Graph of the nonparametric effect ##################

np.graph(fit, which=1, exp=TRUE)

################## Graph of deviance-type residuals ##################

plot(fit)

###################################################################################
####################### gam.data - a Power-exponential model   ####################
###################################################################################

data("gam.data", package="gam")

fit <- ssym.l(y~psp(x),data=gam.data,family="Powerexp",xi=-0.5)
summary(fit)

################## Graph of the nonparametric effect ##################

np.graph(fit, which=1)

################## Graph of deviance-type residuals ##################

plot(fit)

###################################################################################
######### Personal Injury Insurance Data - a Birnbaum-Saunders-t model   ##########
###################################################################################

data("Claims", package="ssym")
fit <- ssym.l(log(total) ~ op_time | op_time, data=Claims,
              family='Sinh-t', xi=c(0.1,4))
summary(fit)

################## Plot of deviance-type residuals ##################

plot(fit)

###################################################################################
######### Body Fat Percentage Data - a Birnbaum-Saunders-t model   ##########
###################################################################################

data("ais", package="sn")
fit <- ssym.l(log(Bfat)~1, data=ais, family='Sinh-t', xi=c(4.5,4))
summary(fit)

id <- sort(ais$Bfat, index=TRUE)$ix
par(mfrow=c(1,2))
hist(ais$Bfat[id],xlim=range(ais$Bfat),ylim=c(0,0.1),prob=TRUE,breaks=15,
     col="light gray",border="dark gray",xlab="",ylab="",main="")
par(new=TRUE)
plot(ais$Bfat[id],exp(fit$lpdf[id])/ais$Bfat[id],xlim=range(ais$Bfat),
     ylim=c(0,0.1),type="l",xlab="",ylab="Density",main="Histogram")
	 
plot(ais$Bfat[id],fit$cdfz[id],xlim=range(ais$Bfat),ylim=c(0,1),type="l",
     xlab="",ylab="",main="")
par(new=TRUE)
plot(ecdf(ais$Bfat[id]),xlim=range(ais$Bfat),ylim=c(0,1),verticals=TRUE,
     do.points=FALSE,col="dark gray",ylab="Probability.",xlab="",main="ECDF")

###################################################################################
################### Boston Housing Data - a log-Slash model   ####################
###################################################################################
	 
#data("Boston", package="MASS")
#fit <- ssym.l(log(medv)~psp(lstat)|psp(lstat),data=Boston,family="Slash",xi=1.7)
#summary(fit)
#plot(fit)
	 
}


