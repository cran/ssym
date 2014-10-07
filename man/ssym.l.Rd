\name{ssym.l}
\alias{ssym.l}
\title{Fitting Semi-parametric Symmetric Regression Models}
\description{
\bold{ssym.l} is used to fit a semi-parametric regression model suitable for 
data set analysis in which the conditional distribution of the response is 
symmetric and continuous. In this setup, both location and dispersion
parameters of the response variable distribution are explicitly modeled
through semi-parametric functions, which may be approximated by natural cubic
splines or P-splines.}
\usage{
ssym.l(formula, family, xi, data, epsilon, maxiter, subset, local.influence)
}
\arguments{
  \item{formula}{a symbolic description of the systematic component of the model to be fitted. See details for further information.}
  \item{family}{a description of the error distribution to be used in the model. Supported families include \emph{Normal}, \emph{Student}, \emph{Contnormal},
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
The argument \emph{formula} comprises three parts (separated by the symbols "~" and "|"), namely: observed response variable, predictor of the
location submodel (having identity link) and predictor of the dispersion submodel (having logarithmic link). A non-parametric effect may be
specified in the predictors, either approximated by a natural cubic spline or a P-spline through the functions \emph{ncs()} or \emph{psp()},
respectively.

The iterative estimation process is based on the Fisher scoring and backfitting algorithms. Because some distributions such as Student-t,
contaminated normal, power exponential, slash and symmetric hyperbolic may be obtained as a scale mixture of normal distributions, the
Expectation-Maximization algorithm is applied in those cases to obtain a more efficient iterative process for the parameter estimation. Furthermore,
because the Sinh-t distribution can be obtained as a scale mixture of Sinh-normal distributions, the Expectation-Maximization algorithm
is also applied in that case to obtain a more efficient iterative process for the parameter estimation. The smoothing parameter(s) is(are)
chosen using the unweighted cross-validation score.

The function \emph{ssym.l()} calculates deviance-type residuals for both submodels as well as local influence measures under case-weight and
response perturbation schemes.
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

\references{Vanegas, L.H. and Paula, G.A. (2014a) A Semiparametric Approach for Joint Modeling of Median and Skewness. TEST (accepted)

			Vanegas, L.H. and Paula, G.A. (2014b) Log-symmetric distributions: statistical properties and parameter estimation. Brazilian Journal of Probability and Statistics (submitted)}

\author{Luis Hernando Vanegas <hvanegasp@gmail.com> and Gilberto A. Paula}

\seealso{
\code{\link{ssym.nl}}
}

\examples{
###################################################################################
########### Fraction of Cell Volume Data - a Power-exponential model  ###########
###################################################################################

data("Ovocytes", package="ssym")
fit <- ssym.l(fraction ~ type + psp(time) | type + psp(time), data=Ovocytes,
                family='Powerexp', xi=-0.65, local.influence=TRUE)
summary(fit)

################## Plot of the non-parametric effects ##################

par(mfrow=c(1,2))
np.graph(fit, which=1)
np.graph(fit, which=2, exp=TRUE)

################## Plot of deviance-type residuals ##################

plot(fit)

################## Plot of local influence measures ##################

ilm <- influence.ssym(fit)

###################################################################################
################## Biaxial Fatigue Data - a Birnbaum-Saunders model   #############
###################################################################################

data("Biaxial", package="ssym")
fit <- ssym.l(log(Life) ~ log(Work), family='Sinh-normal', xi=0.41, data=Biaxial, 
              local.influence=TRUE)
summary(fit)

################## Plot of the fitted model ##################

xl <- "Work per cycle"
rx <-range(Biaxial$Work)
ry <- range(Biaxial$Life)
plot(Biaxial$Work,Biaxial$Life,xlim=rx,ylim=ry,type="p",cex=0.3,lwd=3,ylab="",xlab="")
par(new=TRUE)
ids <- sort(Biaxial$Work, index=TRUE)$ix
plot(Biaxial$Work[ids],exp(fitted(fit)$mu[ids]),xlim=rx,ylim=ry,type="l",ylab="Life",xlab=xl)

################## Plot of deviance-type residuals ##################

plot(fit)

################## Plot of local influence measures ##################

ilm <- influence.ssym(fit)
}
