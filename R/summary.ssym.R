summary.ssym <-
function(object, ...){

familia <- object$family
xxi <- object$xi
if(familia=='Normal'){
cat("\n     Family: ",familia,"\n")}
else{
cat("\n     Family: ",familia," (",xxi,")\n")}
cat("Sample size: ",length(object$deviance.mu),"\n\n")

cat(" Quantile of the Weights\n")
temp <- round(quantile(object$v/object$dg),digits=2)
cmat <- cbind(temp[1], temp[2], temp[3], temp[4], temp[5])
colnames(cmat) <- c("0%", "25%", "50%", "75%", "100%")
rownames(cmat) <- ""
printCoefmat(cmat,digits=3)
temp2 <- qqnorm(qnorm(object$cdfz),plot.it=FALSE)
cat("\n Overall goodness of fit: ",round(mean(abs(sort(temp2$x)-sort(temp2$y))),digits=6),"\n")
cat(" -2*log-likelihood: ",round(-2*sum(object$lpdf),digits=3),"\n")

							  
cat("\n\n ************ Location Model ************\n\n")

TAB		 <- cbind(Estimate <- round(object$coefs.mu,digits=5),
				  StdErr <- round(sqrt(diag(object$vcov.mu)),digits=5),
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
rownames(TAB) <- object$filas
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4)
cat("\n Deviance: ",round(sum(object$deviance.mu),digits=3))


cat("\n\n ************ Dispersion Model ************\n")
if(object$l>0){
        TAB		 <- cbind(Estimate <- round(object$coefs.phi[1:object$l],digits=5),
				  StdErr <- round(object$se.phi[1:object$l],digits=5),
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
        colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
        rownames(TAB) <- object$filas2
  cat("\n Parametric part\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4)		
}
if(object$q>0){
	spar <- object$lambda
	gle <- object$gle
	q <- object$q
  cat("\n Nonparametric part\n")
  cat("                 knots: ",q,"\n")
  cat("   smoothing parameter: ",round(spar,digits=4),"\n")
  cat("    degrees of freedom: ",round(gle,digits=2),"\n")
}
cat("\n Deviance: ",round(sum(object$deviance.phi),digits=3),"\n")
}
