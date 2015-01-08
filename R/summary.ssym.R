summary.ssym <-
function(object, ...){

familia <- object$family
xxi <- object$xi
if(familia=='Normal'){
cat("\n     Family: ",familia,"\n")}
else{
    if(familia=='Contnormal' | familia=='Sinh-t')
      cat("\n     Family: ",familia,"(",xxi[1],",",xxi[2],")\n")
	else cat("\n     Family: ",familia,"(",xxi[1],")\n")
}

cat("Sample size: ",length(object$deviance.mu),"\n")

cat(" Quantile of the Weights\n")
temp <- round(quantile(object$weights),digits=2)
cmat <- cbind(temp[1], temp[2], temp[3], temp[4], temp[5])
colnames(cmat) <- c("0%", "25%", "50%", "75%", "100%")
rownames(cmat) <- ""
printCoefmat(cmat,digits=3)
gle <- 0
							  
cat("\n ************************** Median/Location submodel **************************\n")
if(object$p>0){
gle <- gle + object$p
TAB		 <- cbind(Estimate <- round(object$coefs.mu[1:object$p],digits=5),
				  StdErr <- round(sqrt(diag(object$vcov.mu))[1:object$p],digits=5),
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
rownames(TAB) <- object$filas
  cat(" **** Parametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4, signif.legend=FALSE)}
if(object$qm>0){
	spar <- object$lambda.mu
	gle <- gle + object$dfe.mu
	q <- object$qm
	B <- cbind(1,-diag(q-1))
	f <- object$coefs.mu[(object$p+1):(object$p+q)]
	v <- (B%*%object$vcov.mu[(object$p+1):(object$p+q),(object$p+1):(object$p+q)]%*%t(B))
	v2 <- eigen(v,symmetric=TRUE)
  cat("\n **** Nonparametric component\n")
  cat("                 knots: ",q,"\n")
  cat("    degrees of freedom: ",round(object$dfe.mu,digits=3),"\n")
  if(min(v2$values)/max(v2$values) > 1e-15)
  cat("               p-value: ",format.pval(1-pchisq(t(B%*%f)%*%solve(v)%*%(B%*%f),q-1),digits=3),"\n")
}
cat("\n **** Deviance: ",round(sum(object$deviance.mu),digits=2),"\n")
gle2 <- 0
cat(" *********************** Skewness/Dispersion submodel ***********************\n")
if(object$l>0){
gle2 <- gle2 + object$l
        TAB		 <- cbind(Estimate <- round(object$coefs.phi[1:object$l],digits=5),
				  StdErr <- round(sqrt(diag(object$vcov.phi))[1:object$l],digits=5),
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
        colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
        rownames(TAB) <- object$filas2
  cat(" **** Parametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4, signif.legend=FALSE)		
}
if(object$q>0){
	spar <- object$lambda.phi
	gle2 <- gle2 + object$dfe.phi
	q <- object$q
	B <- cbind(1,-diag(q-1))
	f <- object$coefs.phi[(object$l+1):(object$l+q)]
	v <- (B%*%object$vcov.phi[(object$l+1):(object$l+q),(object$l+1):(object$l+q)]%*%t(B))
	v2 <- eigen(v,symmetric=TRUE)
  cat("\n **** Nonparametric component\n")
  cat("                 knots: ",q,"\n")
  cat("    degrees of freedom: ",round(object$dfe.phi,digits=3),"\n")
  if(!is.complex(v2$values)){
	  if(min(v2$values)/max(v2$values) > 1e-15)  
	  cat("               p-value: ",format.pval(1-pchisq(t(B%*%f)%*%solve(v)%*%(B%*%f),q-1),digits=3),"\n")}
}
cat("\n **** Deviance: ",round(sum(object$deviance.phi),digits=2),"\n")
cat(" ****************************************************************************\n")
temp2 <- qqnorm(qnorm(object$cdfz),plot.it=FALSE)
cat(" Overall goodness-of-fit statistic: ",round(mean(abs(sort(temp2$x)-sort(temp2$y))),digits=6),"\n")
cat("                 -2*log-likelihood: ",round(-2*sum(object$lpdf),digits=3),"\n")
cat("                               AIC: ",round(-2*sum(object$lpdf) + 2*(gle+gle2),digits=3),"\n")
cat("                               BIC: ",round(-2*sum(object$lpdf) + log(length(object$lpdf))*(gle+gle2),digits=3),"\n")
}
