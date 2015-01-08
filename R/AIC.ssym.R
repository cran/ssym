AIC.ssym <-
function(object, ...){
	gle <- 0
	if(object$p>0) gle <- gle + object$p
	if(object$qm>0)	gle <- gle + object$dfe.mu
	if(object$l>0) gle <- gle + object$l
	if(object$q>0) gle <- gle + object$dfe.phi

    AIC <- round(-2*sum(object$lpdf) + 2*(gle), digits=3)
    y <- object$z.hat*sqrt(object$phi.fitted) + object$mu.fitted
    attr(AIC,"log") <- round(-2*sum(object$lpdf) + 2*(gle) + 2*sum(y), digits=3)
    AIC
}
