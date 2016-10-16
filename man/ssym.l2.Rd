\name{ssym.l2}
\alias{ssym.l2}
\title{Fitting Censored Semi-parametric Log-symmetric Regression Models}
\description{
\bold{ssym.l2} is used to fit a semi-parametric regression model suitable for 
analysis of data sets in which the response variable is continuous, strictly
positive, asymmetric and right-censored. Under this setup, both median and skewness of the response
variable distribution are explicitly modeled through semi-parametric functions,
whose nonparametric components may be approximated by natural cubic splines or P-splines.}
\usage{
ssym.l2(formula, family, xi, data, epsilon, maxiter, subset, link.mu, link.phi,
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
  \item{subset}{an optional expression specifying a subset of individuals to be used in the fitting process.}
  \item{link.mu}{an optional character that specifies the link function of the median submodel.}
  \item{link.phi}{an optional character that specifies the link function of the skewness submodel.}
  \item{local.influence}{logical. If TRUE, local influence measures under two perturbation schemes are calculated. Default is \code{FALSE}.}
  \item{spec}{character. The smoothing parameter is estimated by minimizing a overall goodness-of-fit criterion such as AIC or BIC. \code{spec} is an optional
              string to specify the goodness-of-fit measure to be used. Default value is \code{AIC}.}
  \item{std.out}{logical. If FALSE, just a reduced set of attributes is returned by the model-fitting function. Default is \code{TRUE}}.
}
\details{
The argument \emph{formula} comprises of three parts (separated by the symbols "~" and "|"), namely: event status and observed response variable (in log-scale) in a object
of class \emph{Surv}, predictor of the median submodel (having logarithmic link) and predictor of the skewness (or the relative dispersion) submodel (having logarithmic link).
An arbitrary number of nonparametric effects may be specified in the predictors. These effects are specified to be approximated by natural cubic splines or P-splines using the functions
\code{ncs()} or \code{psp()}, respectively.

The iterative estimation process is based on the Gauss-Seidel, Newton-Raphson and backfitting algorithms. The
smoothing parameter is chosen by minimizing the AIC or BIC criteria.

The function \code{ssym.l2()} calculates overall goodness-of-fit statistics, deviance-type residuals for both submodels, as well as local influence measures
under the case-weight and response perturbation schemes.
}

\value{
\item{theta.mu}{a vector of parameter estimates associated with the median submodel.}
\item{theta.phi}{a vector of parameter estimates associated with the skewness (or the relative dispersion) submodel.}
\item{vcov.mu}{approximate variance-covariance matrix associated with the median submodel.}
\item{vcov.phi}{approximate variance-covariance matrix associated with the skewness (or the relative dispersion) submodel.}
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
\code{\link{ssym.nl}, \link{ssym.l}}
}


\examples{

###################################################################################
################ Lung Cancer Trial - a log-Student model ##########################
###################################################################################
#data("veteran", package="survival")
#fit <- ssym.l2(Surv(log(time), status) ~ karno| karno, data = veteran,
#              family="Student", xi=4.5)
#summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,3,10)
#
################## Graph of deviance-type residuals ##################
#plot(fit)

####################################################################################
########## Primary biliary cirrhosis - a Power-exponential model ###################
####################################################################################
# data("pbc", package="survival")
# pbc2 <- data.frame(pbc[!is.na(pbc$edema) & !is.na(pbc$stage) & !is.na(pbc$bili),])
#
# fit <- ssym.l2(Surv(log(time),ifelse(status>=1,1,0) ) ~ factor(edema) +
#                stage + ncs(bili), data = pbc2, family="Powerexp",
#                xi=0.47, local.influence=TRUE)
# summary(fit)
#
########################### Extra parameter ###########################
#extra.parameter(fit,c(0.6,3),c(0.9,5))
#
################## Graph of the nonparametric effect ##################
#np.graph(fit, which=1, exp=TRUE)
#
################## Graph of deviance-type residuals ##################
#plot(fit)

####################################################################################
################### Myeloma - a Birnbaum-Saunders model ###########################
####################################################################################
# data("myeloma", package="ssym")
#
#fit <- ssym.l2(Surv(log(t),1-event) ~ x1 + x2 + x5| -1 + x3, data=myeloma,
#       family="Sinh-normal", xi=1.8)
#summary(fit)
#
################## Graph of deviance-type residuals ##################
#plot(fit)

####################################################################################
################ Baboons Data - a log-power-exponential model   ####################
####################################################################################
########################## left-censored observations ##############################
####################################################################################

#data("Baboons", package="ssym")
#fit <- ssym.l2(Surv(-log(t),1-cs) ~ 1, data=Baboons, family="Powerexp", xi=-0.35)
#
#summary(fit)
################## Graph of deviance-type residuals ##################
#plot(fit)

}