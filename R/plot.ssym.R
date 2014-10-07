plot.ssym <-
function(x, ...){
par(mfrow=c(2,2))
plot(x$z.hat, x$weights, main = "Individual-specific Weights", cex=0.3, lwd=3, xlab="Ordinary residuals", ylab="Weight")
qqnorm(qnorm(x$cdfz), main="Overall goodness-of-fit statistic", cex=0.3, lwd=3, xlab="Quantiles of N(0,1)", ylab="Overall residuals")
abline(0,1)
res.dev.mu <- sqrt(x$deviance.mu)*ifelse(x$z.hat>=0,1,-1)   
ry <- c(min(res.dev.mu,-3.5),max(res.dev.mu,3.5))               
plot(x$mu.fitted,res.dev.mu, ylim=ry, cex=0.3, lwd=3, main="Median/location submodel", ylab="Deviance-type residual", xlab="Fitted values")                                 
abline(h=-3,lty=3)                                              
abline(h=+3,lty=3)                                              
res.dev.phi <- sqrt(x$deviance.phi)*ifelse(x$z.hat>=0,1,-1) 
ry <- c(min(res.dev.phi,-3.5),max(res.dev.phi,3.5))             
plot(x$phi.fitted,res.dev.phi, ylim=ry, cex=0.3, lwd=3, main="Skewness/dispersion submodel", ylab="Deviance-type residual", xlab="Fitted values")                                
abline(h=-3,lty=3)                                              
abline(h=+3,lty=3)
}
