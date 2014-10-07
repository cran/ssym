print.ssym <-
function(x, ...){

familia <- x$family
xxi <- x$xi
if(familia=='Normal'){
cat("\n     Family: ",familia,"\n")}
else{
    if(familia=='Contaminated Normal' | familia=='Sinh-t')
      cat("\n     Family: ",familia,"(",xxi[1],",",xxi[2],")\n")
	else cat("\n     Family: ",familia,"(",xxi[1],")\n")
}

cat("\n ************************** Median/Location submodel **************************\n")
if(x$p>0){
TAB		 <- cbind(Estimate <- round(x$coefs.mu[1:x$p],digits=5),
				  StdErr <- round(sqrt(diag(x$vcov.mu))[1:x$p],digits=5),
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
rownames(TAB) <- x$filas
  cat(" **** Parametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4, signif.legend=FALSE)}
if(x$qm>0){
	spar <- x$lambda.mu
	q <- x$qm
  cat("\n **** Nonparametric component\n")
  cat("                 knots: ",q,"\n")
  cat("   smoothing parameter: ",round(spar,digits=3),"\n")
}

cat("\n\n *********************** Skewness/Dispersion submodel ***********************\n")
if(x$l>0){
        TAB		 <- cbind(Estimate <- round(x$coefs.phi[1:x$l],digits=5),
				  StdErr <- round(sqrt(diag(x$vcov.phi))[1:x$l],digits=5),
				  tval <- Estimate/StdErr,
				  p.value <- 2*pnorm(-abs(tval)))
        colnames(TAB) <- c("Estimate", "Std.Err", "z-value", "Pr(>|z|)")
        rownames(TAB) <- x$filas2
  cat(" **** Parametric component\n\n")		
printCoefmat(TAB, P.values=TRUE, has.Pvalue=TRUE,digits=4, signif.legend=FALSE)		
}
if(x$q>0){
	spar <- x$lambda.phi
	q <- x$q
  cat("\n **** Nonparametric component\n")
  cat("                 knots: ",q,"\n")
  cat("   smoothing parameter: ",round(spar,digits=3),"\n")
}
}
