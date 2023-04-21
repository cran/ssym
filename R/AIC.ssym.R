#' \name{AIC.ssym}
#' \alias{AIC.ssym}
#' \title{AIC.ssym}
#' \description{
#'   \bold{AIC.ssym} calculates the goodness-of-fit statistic AIC from an object of class ``"ssym".}
#'
AIC.ssym <-
function(object, ...){
	gle <- sum(object$gle.mu) + sum(object$gle.phi)

    AIC <- round(-2*sum(object$lpdf) + 2*(gle), digits=3)
    y <- object$z_es*sqrt(object$phi.fitted) + object$mu.fitted
	if(object$censored==FALSE) attr(AIC,"log") <- round(-2*sum(object$lpdf) + 2*(gle) + 2*sum(y), digits=3)
	else attr(AIC,"log") <- round(-2*sum(object$lpdf) + 2*(gle) + 2*sum(y[object$event==0]), digits=3)
    AIC
}
