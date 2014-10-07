residuals.ssym <-
function(object, ...){
list(mu=sqrt(object$deviance.mu)*ifelse(object$z.hat>=0,1,-1), phi=sqrt(object$deviance.phi)*ifelse(object$z.hat>=0,1,-1), overall=qnorm(object$cdfz), ordinary=object$z.hat)}
