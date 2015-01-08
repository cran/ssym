BIC.ssym <-
function(object, ...){
	gle <- 0
	if(object$p>0) gle <- gle + object$p
	if(object$qm>0)	gle <- gle + object$dfe.mu
	if(object$l>0) gle <- gle + object$l
	if(object$q>0) gle <- gle + object$dfe.phi

    BIC <- round(-2*sum(object$lpdf) + log(length(object$mu.fitted))*(gle), digits=3)
    y <- object$z.hat*sqrt(object$phi.fitted) + object$mu.fitted
    attr(BIC,"log") <- round(-2*sum(object$lpdf) + log(length(object$mu.fitted))*(gle) + 2*sum(y), digits=3)
    BIC
}
