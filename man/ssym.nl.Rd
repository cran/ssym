\name{ssym.nl}
\alias{ssym.nl}
\title{Fitting Semi-parametric Log-symmetric Regression Models}
\description{
\bold{ssym.nl} is used to fit a semi-parametric regression model suitable for 
analysis of data sets in which the response variable is continuous, strictly
positive, and asymmetric. Under this setup, both median and skewness
of the response variable distribution are explicitly modeled,
the median using a nonlinear function and the skewness through semi-parametric
functions, whose nonparametric components may be approximated by natural cubic splines or P-splines.}
\usage{
ssym.nl(formula, start, family, xi, data, epsilon, maxiter, subset, link.phi,
       local.influence, spec, std.out)
}
\arguments{
  \item{formula}{a symbolic description of the systematic component of the model to be fitted. See details for further information.}
  \item{start}{a named numeric vector of starting estimates for the parameters in the specified nonlinear function.}  
  \item{family}{a description of the (log) error distribution to be used in the model. Supported families include \code{Normal}, \code{Student}, \code{Contnormal},
               \code{Powerexp}, \code{Hyperbolic}, \code{Slash}, \code{Sinh-normal} and \code{Sinh-t}, which correspond to normal, Student-t, contaminated normal,
			    power exponential, symmetric hyperbolic, slash, sinh-normal and sinh-t distributions, respectively.}
  \item{xi}{a numeric value or numeric vector that represents the extra parameter of the specified error distribution.}
  \item{data}{an optional data frame, list or environment containing the variables in the model.}
  \item{epsilon}{an optional positive value, which represents the convergence criterion. Default value is 1e-07.}
  \item{maxiter}{an optional positive integer giving the maximal number of iterations for the estimating process. Default value is 1e03.}
  \item{subset}{an optional expression that specifies a subset of individuals to be used in the fitting process.}
  \item{link.phi}{an optional character that specifies the link function of the skewness submodel.}
  \item{local.influence}{logical. If TRUE, local influence measures under two perturbation schemes are calculated. Default is \code{FALSE}.}
  \item{spec}{character. The smoothing parameter is estimated by minimizing a overall goodness-of-fit criterion such as AIC or BIC. \code{spec} is an optional
              string to specify the goodness-of-fit measure to be used. Default value is \code{AIC}.}
  \item{std.out}{logical. If FALSE, just a reduced set of attributes is returned by the model-fitting function. Default is \code{TRUE}}.
}
\details{
The argument \emph{formula} comprises of three parts (separated by the symbols "~" and "|"), namely: observed response variable in log-scale, predictor of the
median submodel (having logarithmic link) and predictor of the skewness (or the relative dispersion) submodel (having logarithmic link). An arbitrary number of
nonparametric effects may be specified in the predictor of the skewness submodel. These effects are specified to be approximated by natural cubic splines or P-splines using the functions
\code{ncs()} or \code{psp()}, respectively.

The iterative estimation process is based on the Fisher scoring and backfitting algorithms. Because some distributions such as log-Student-t,
log-contaminated-normal, log-power-exponential, log-slash and log-hyperbolic may be obtained as a power mixture of the log-normal distribution, the
expectation-maximization (EM) algorithm is applied in those cases to obtain a more efficient iterative process of parameter estimation. Furthermore,
because the Birnbaum-Saunders-t distribution can be obtained as a scale mixture of the Birnbaum-Saunders distribution, the expectation-maximization algorithm
is also applied in this case to obtain a more efficient iterative process of parameter estimation. The smoothing parameter is chosen by minimizing
the AIC or BIC criteria.

The function \code{ssym.nl()} calculates overall goodness-of-fit statistics, deviance-type residuals for both submodels, as well as local influence measures
under the case-weight and response perturbation schemes.
}

\value{
\item{theta.mu}{a vector of parameter estimates associated with the median submodel.}
\item{theta.phi}{a vector of parameter estimates associated with the skewness (or the relative dispersion) submodel.}
\item{vcov.mu}{approximate variance-covariance matrix associated with the median submodel.}
\item{vcov.phi}{approximate variance-covariance matrix associated with the skewness (or the relative dispersion) submodel.}
\item{weights}{final weights of the iterative process.}
\item{lambdas.phi}{estimate of the smoothing parameter(s) associated with the nonparametric part of the skewness (or the relative dispersion) submodel.}
\item{gle.mu}{degrees of freedom associated with the nonparametric part of the median submodel.}
\item{gle.phi}{degrees of freedom associated with the nonparametric part of the skewness (or the relative dispersion) submodel.}
\item{deviance.mu}{a vector with the individual contributions to the deviance associated with the median submodel.}
\item{deviance.phi}{a vector with the individual contributions to the deviance associated with the skewness (or the relative dispersion) submodel.}
\item{mu.fitted}{a vector with the fitted values of the (in log-scale) median submodel.}
\item{phi.fitted}{a vector with the fitted values of the skewness (or the relative dispersion) submodel.}
\item{lpdf}{a vector of individual contributions to the log-likelihood function.}
}

\references{Vanegas, L.H. and Paula, G.A. (2015a) A Semiparametric Approach for Joint Modeling of Median and Skewness. TEST 24, 110-135.

			Vanegas, L.H. and Paula, G.A. (2015b) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics (to appear)}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\seealso{
\code{\link{ssym.l}, \link{ssym.l2}}
}

\examples{

###################################################################################
######### Ultrasonic Calibration Data - a log-contaminated-normal model ###########
###################################################################################
#data("Chwirut1", package="NISTnls")
#fit <- ssym.nl(log(y) ~ -b1*x-log(b2 + b3*x) | x, start=c(b1=0.15,b2=0.005,b3=0.012),
#        data=Chwirut1, family='Contnormal', xi=c(0.68,0.1), local.influence=TRUE)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,c(0.4,0.08),c(0.9,0.11))
#
################## Graph of deviance-type residuals ##################
#plot(fit)
#
################## Simulated envelopes ##################
#envelope(fit)
#
################## Graph of local influence measures ##################
#out <- influence.ssym(fit)

###################################################################################
################## Biaxial Fatigue Data - a Birnbaum-Saunders model   #############
###################################################################################
#data("Biaxial", package="ssym")
#fit <- ssym.nl(log(Life) ~ b1*Work^b2, start=c(b1=16, b2=-0.25),
#                data=Biaxial, family='Sinh-normal', xi=1.54)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,1.3,1.8)

################## Graph of deviance-type residuals ##################
#plot(fit)

################## Simulated envelopes ##################
#envelope(fit)

###################################################################################
################## European rabbits Data - a log-normal model   #############
###################################################################################
#data("Erabbits", package="ssym")
#fit <- ssym.nl(log(wlens) ~ b1 - b2/(b3 + age) | psp(age), start=c(b1=5,
#               b2=130, b3=36), data=Erabbits, family='Normal')
#summary(fit)
#
################## Graph of the nonparametric effect ##################
#np.graph(fit, which=2, exp=TRUE)
#
################## Graph of deviance-type residuals ##################
#plot(fit)

###################################################################################
################### Metsulfuron Data - a log-slah model ########################
###################################################################################
#data("M4", package="nlreg")
#fit <- ssym.nl(log(area) ~ log(b1+(b2-b1)/(1+(dose/b4)^b3)) | dose, data=M4,
#       start = c(b1 = 138, b2 = 2470, b3 = 2, b4 = 0.07), family="Slash", xi=1.4)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,1.0,1.6)
#
################## Graph of deviance-type residuals ##################
#plot(fit)
#
################## Simulated envelopes ##################
#envelope(fit)

###################################################################################
################# Blood flow Data - a log-power-exponential model #################
###################################################################################
#data("la", package="gamlss.nl")
#fit <- ssym.nl(log(PET60) ~ log(bflow) + log(1+b1*exp(-b2/bflow)) | bflow,
#       data=la, start=c(b1=-0.6,b2=98), family="Powerexp", xi=-0.45)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,-0.5,0)

################## Graph of deviance-type residuals ##################
#plot(fit)

###################################################################################
######### Gross Domestic Product per capita Data - a Birnbaum-Saunders model ######
###################################################################################
#data("gdp", package="ssym")
#fit <- ssym.nl(log(gdp2010) ~ b1, start=c(b1=mean(log(gdp$gdp2010))), data=gdp, 
#               family='Sinh-normal', xi=2.2)
#summary(fit)
########################### Extra parameter ###########################
#extra.parameter(fit,0.5,3)

################## Plot of the fitted model ##################
#id <- sort(gdp$gdp2010, index=TRUE)$ix
#par(mfrow=c(1,2))
#hist(gdp$gdp2010[id],xlim=range(gdp$gdp2010),ylim=c(0,0.00025),prob=TRUE,
#     breaks=200,col="light gray",border="dark gray",xlab="",ylab="",main="")
#par(new=TRUE)
#plot(gdp$gdp2010[id],exp(fit$lpdf[id])/gdp$gdp2010[id],xlim=range(gdp$gdp2010),
#     ylim=c(0,0.00025),type="l",xlab="",ylab="Density",main="Histogram")
#	 
#plot(gdp$gdp2010[id],fit$cdfz[id],xlim=range(gdp$gdp2010),ylim=c(0,1),type="l",
#     xlab="",ylab="",main="")
#par(new=TRUE)
#plot(ecdf(gdp$gdp2010[id]),xlim=range(gdp$gdp2010),ylim=c(0,1),verticals=TRUE,
#     do.points=FALSE,col="dark gray",ylab="Probability.",xlab="",main="ECDF")

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
#
################## Graph of deviance-type residuals ##################
#plot(fit)

###################################################################################
################ Daphnia Data - a log-power-exponential model #####################
###################################################################################
#data("daphnia", package="nlreg")
#fit <- ssym.nl(log(time) ~ log(b1+(b2-b1)/(1+(conc/b4)^b3)) | ncs(conc), data=daphnia,
#       start = c(b1=0, b2=50 , b3=2, b4=0.2), family="Powerexp", xi=-0.4)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,-0.6,0)
#
################## Graph of deviance-type residuals ##################
#plot(fit)
#
################## Simulated envelopes ##################
#envelope(fit)

}