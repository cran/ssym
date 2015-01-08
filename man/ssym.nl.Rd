\name{ssym.nl}
\alias{ssym.nl}
\title{Fitting Semi-parametric Log-symmetric Regression Models}
\description{
\bold{ssym.nl} is used to fit a semi-parametric regression model suitable for 
analysis of data sets in which the response variable is continuous, strictly
positive, and asymmetric. In this setup, both median and skewness
of the response variable distribution are explicitly modeled,
the median using a nonlinear function and the skewness using a semi-parametric
function, which may be approximated by a natural cubic spline or a P-spline.}
\usage{
ssym.nl(formula, start, family, xi, data, epsilon, maxiter, subset, local.influence)
}
\arguments{
  \item{formula}{a symbolic description of the systematic component of the model to be fitted. See details for further information.}
  \item{start}{a named numeric vector of starting estimates for the parameters in the specified nonlinear function.}  
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

The function \emph{ssym.nl()} calculates deviance-type residuals for both submodels as well as local influence measures under case-weight and
response perturbation schemes.
}

\value{
\item{coefs.mu}{a vector of parameter estimates associated with the median submodel.}
\item{coefs.phi}{a vector of parameter estimates associated with the skewness submodel.}
\item{vcov.mu}{approximate variance-covariance matrix associated with the median submodel.}
\item{vcov.phi}{approximate variance-covariance matrix associated with the skewness submodel.}
\item{weights}{final weights of the iterative process.}
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
\code{\link{ssym.l}}
}

\examples{

###################################################################################
######### Ultrasonic Calibration Data - a log-contaminated-normal model ###########
###################################################################################

data("Chwirut1", package="NISTnls")
sv <- c(b1=0.15, b2=0.005, b3=0.012)
fit <- ssym.nl(log(y) ~ -b1*x - log(b2 + b3*x) | psp(x), start=sv,
       data=Chwirut1, family='Contnormal', xi=c(0.6,0.1), local.influence=TRUE)
summary(fit)

################## Graph of the nonparametric effect ##################

np.graph(fit, which=2, exp=TRUE)

################## Graph of deviance-type residuals ##################

plot(fit)

################## Graph of local influence measures ##################

ilm <- influence.ssym(fit)

###################################################################################
################## Biaxial Fatigue Data - a Birnbaum-Saunders model   #############
###################################################################################

data("Biaxial", package="ssym")
sv <- c(b1=16, b2=-0.25)
fit <- ssym.nl(log(Life) ~ b1*Work^b2, start=sv, data=Biaxial,
               family='Sinh-normal', xi=1.54)
summary(fit)

################## Graph of deviance-type residuals ##################

plot(fit)

###################################################################################
################## European rabbits Data - a log-normal model   #############
###################################################################################

data("Erabbits", package="ssym")
fit <- ssym.nl(log(wlens) ~ b1 - b2/(b3 + age) | psp(age), start=c(b1=5.6,
               b2=128, b3=36.4), data=Erabbits, family='Normal')
summary(fit)

################## Graph of the nonparametric effect ##################

np.graph(fit, which=2, exp=TRUE)

################## Graph of deviance-type residuals ##################

plot(fit)

###################################################################################
######### Gross Domestic Product per capita Data - a Birnbaum-Saunders model ######
###################################################################################

data("gdp", package="ssym")
fit <- ssym.nl(log(gdp2010) ~ b1, start=c(b1=mean(log(gdp$gdp2010))), data=gdp, 
               family='Sinh-normal', xi=2.2)
summary(fit)

################## Plot of the fitted model ##################

id <- sort(gdp$gdp2010, index=TRUE)$ix
par(mfrow=c(1,2))
hist(gdp$gdp2010[id],xlim=range(gdp$gdp2010),ylim=c(0,0.00025),prob=TRUE,
     breaks=200,col="light gray",border="dark gray",xlab="",ylab="",main="")
par(new=TRUE)
plot(gdp$gdp2010[id],exp(fit$lpdf[id])/gdp$gdp2010[id],xlim=range(gdp$gdp2010),
     ylim=c(0,0.00025),type="l",xlab="",ylab="Density",main="Histogram")
	 
plot(gdp$gdp2010[id],fit$cdfz[id],xlim=range(gdp$gdp2010),ylim=c(0,1),type="l",
     xlab="",ylab="",main="")
par(new=TRUE)
plot(ecdf(gdp$gdp2010[id]),xlim=range(gdp$gdp2010),ylim=c(0,1),verticals=TRUE,
     do.points=FALSE,col="dark gray",ylab="Probability.",xlab="",main="ECDF")


###################################################################################
################# Blood flow Data - a log-power-exponential model #################
###################################################################################

#data("la", package="gamlss.nl")
#fit <- ssym.nl(log(PET60) ~ log(bflow) + log(1+b1*exp(-b2/bflow)) | bflow,
#               data=la, start=c(b1=-0.6,b2=98), family="Powerexp", xi=-0.45)
#summary(fit)
#plot(fit)


###################################################################################
############# Australian Institute of Sport Data - a log-normal model #############
###################################################################################

#data("ais", package="sn")
#sex <- ifelse(ais$sex=="male",1,0)
#ais2 <- data.frame(BMI=ais$BMI,LBM=ais$LBM,sex)
#start = c(b1=7, b2=0.3, b3=2, b4=0)
#fit <- ssym.nl(log(BMI) ~ log(b1 + b2*LBM + b3*sex + b4*LBM*sex) | sex + LBM,
#               data=ais2, start=start, family="Normal")
#summary(fit)
#plot(fit)

}
