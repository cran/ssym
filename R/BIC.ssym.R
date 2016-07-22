BIC.ssym <-
function(object, ...){
	gle <- sum(object$gle.mu) + sum(object$gle.phi)

    BIC <- round(-2*sum(object$lpdf) + log(length(object$mu.fitted))*(gle), digits=3)
    y <- object$z_es*sqrt(object$phi.fitted) + object$mu.fitted
	if(object$censored==FALSE) attr(BIC,"log") <- round(-2*sum(object$lpdf) + log(length(y))*(gle) + 2*sum(y), digits=3)
	else attr(BIC,"log") <- round(-2*sum(object$lpdf) + log(length(y))*(gle) + 2*sum(y[object$event==0]), digits=3)
    BIC
}
