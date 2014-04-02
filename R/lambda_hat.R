lambda_hat <-
function(response,x,xx,lambda,plot){
	z <- splinek(x,xx)
	N <- z$N
	K <- z$K
	fr <- function(eta){
		  lambda <- exp(eta)
	      H <- N%*%solve(t(N)%*%N + lambda*K)%*%t(N)
	      h <- diag(H)
	      ajust <- H%*%response
	      CV<- mean(((response-ajust)/(1-h))^2)
	CV
	}
	start <- lambda
	salida <- nlminb(log(start),fr)
	if(!missingArg(plot)){
	  lambdas <- seq(salida$par*0.9,salida$par*1.1,length=50)
	  cves <- matrix(0,length(lambdas),1)
	  for(i in 1:length(lambdas)) cves[i] <- fr(lambdas[i])
	  plot(exp(lambdas),cves,type="l",xlab="smoothing parameter",ylab="cross-validation score",main="Behavior of the cross-validation score")
	  abline(v=exp(salida$par),col="blue",lty=3)
	}
	list(value=salida$value, lambda_hat=exp(salida$par))
}
