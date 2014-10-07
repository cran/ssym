lambda.hat <-
function(response,xx,lambda,type,plot){
	if(type!=1 & type!=2) stop("The type value must be 1 for ncs or 2 for psp!!",call.=FALSE)
	if(type==1) z <- ncs(xx)
	if(type==2) z <- psp(xx)
	N <- attr(z,"N")
	K <- attr(z,"K")
	fr <- function(eta){
		  lambda <- exp(eta)
	      H <- N%*%solve(t(N)%*%N + lambda*K)%*%t(N)
	      h <- diag(H)
	      ajust <- H%*%(response)
	      CV<- mean((((response-ajust)/(1-h))^2))
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
	list(value=fr(salida$par), lambda_hat=exp(salida$par))
}
