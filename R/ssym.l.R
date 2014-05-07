ssym.l <-
function(response, formula.mu, ncs.mu, start.lambda.mu, lambda.mu, formula.phi, ncs.phi, start.lambda.phi, lambda.phi, family, xi, epsilon, maxiter, subset, local.influence){

if(family!="Normal" & family!="Slash" & family!="Hyperbolic" & family!="Sinh-t" & family!="Sinh-normal" & family!="Contnormal" & family!="Powerexp" & family!="Student")
stop("family of distributions specified by the user is not supported!!",call.=FALSE)

if(family=="Slash" | family=="Hyperbolic" | family=="Sinh-t" | family=="Sinh-normal" | family=="Contnormal" | family=="Powerexp" | family=="Student"){
  if(missingArg(xi)) stop("for the family of distributions specified by the user an extra parameter is required!!", call.=FALSE) 
}

if(missingArg(epsilon)) epsilon <- 0.0000001
if(missingArg(maxiter)) maxiter <- 1000
if(missingArg(start.lambda.phi)) start.lambda.phi <- 1
if(missingArg(start.lambda.mu)) start.lambda.mu <- 1
if(missingArg(subset)) subset <- 1:length(response)
y <- response[subset]
n <- length(y)

if(missingArg(ncs.mu)){
  qm <- 0
  if(missingArg(formula.mu)){
    p <- 1
    X <- matrix(1,length(y),1)
	colnames(X) <- "Intercept"
  }
  else{
  	mf <- model.frame(formula=formula.mu)
  	X <- as.matrix(model.matrix(attr(mf, "terms"), data=mf))[subset,]
	if(!is.matrix(X)){
	   XX <- matrix(0,n,1)
	   XX[,1] <- X
	   X <- XX
	}
  	p <- ncol(X)
  }
}
else{
   xxm <- ncs.mu[subset]
   xm <- as.numeric(levels(factor(xxm)))
   qm <- length(xm)
   if(missingArg(formula.mu)){
    p <- 0
   }
   else{
	mf <- model.frame(formula=formula.mu)
	WW <- model.matrix(attr(mf, "terms"), data=mf)
	X <- as.matrix(WW[subset,2:ncol(WW)])
	colnames(X) <- colnames(model.matrix(attr(mf, "terms"), data=mf))[2:ncol(WW)]
	p <- ncol(X)
   }
}
if(p>0){
	mu <- function(vP){X%*%vP[1:p]}
	GradD <- function(vP){X}
	start <- solve(t(X)%*%X)%*%t(X)%*%y
	beta0 <- start
}else{
	start <- rep(mean(y),qm)
}



if(missingArg(ncs.phi)){
  q <- 0
  if(missingArg(formula.phi)){
    l <- 1
    W <- matrix(1,length(y),1)
	colnames(W) <- "Intercept"
  }
  else{
  	mf <- model.frame(formula=formula.phi)
  	W <- as.matrix(model.matrix(attr(mf, "terms"), data=mf))[subset,]
	if(!is.matrix(W)){
	   WW <- matrix(0,n,1)
	   WW[,1] <- W
	   W <- WW
	}
	l <- ncol(W)
  }
}
else{
   xx <- ncs.phi[subset]
   x <- as.numeric(levels(factor(xx)))
   q <- length(x)
   if(missingArg(formula.phi)){
    l <- 0
   }
   else{
	mf <- model.frame(formula=formula.phi)
	WW <- model.matrix(attr(mf, "terms"), data=mf)
	W <- as.matrix(WW[subset,2:ncol(WW)])
	l <- ncol(W)
   }
}
subset <- 1:n


if(family=="Normal"){
	xi <- 0
	v <- function(z){
	     rep(1,length(z))
	}
	vp <- function(z){
	      rep(0,length(z))
	}
	dg <- 1
    fg <- 3
	deviance.mu <- function(z){z^2}
	deviance.phi <- function(z){z^2-1-log(z^2)}
	cdfz <- function(z){pnorm(z)}
	lpdf <- function(z,phi){
	       log(exp(-z^2/2)/sqrt(2*pi)) -(1/2)*log(phi)
    }
	xix <- 1
}
if(family=="Student"){
	if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	family="Student-t"
	nu <- xi[1]
	v <- function(z){
		 (nu + 1)/(nu + z^2)
	}
	vp <- function(z){
		 -2*z*(nu + 1)/((nu + z^2)^2)
	}
	dg <- (nu + 1)/(nu + 3)
    fg <- 3*(nu + 1)/(nu + 3)
	deviance.mu <- function(z){(nu+1)*log(1 + z^2/nu)}
	deviance.phi <- function(z){(nu+1)*log((nu + z^2)/(nu + 1)) -log(z^2)}
	cdfz <- function(z){pt(z,nu)}
	lpdf <- function(z,phi){
	       log((gamma((nu+1)/2)/(sqrt(nu*pi)*gamma(nu/2)))*(1 + z^2/nu)^(-(nu+1)/2)) -(1/2)*log(phi)

    }
	xix <- nu/(nu-2)
	if(nu<=2) xix <- as.null(xix)
}
if(family=="Contnormal"){
	if(xi[1]<=0  | xi[1]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
	if(xi[2]<=0  | xi[2]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
	family="Contaminated Normal"	
	eta <- xi[1:2]
	v <- function(z){
		 (eta[2]^(3/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))/(eta[2]^(1/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))
		 
	}
	vp <- function(z){
		  eta[2]^(1/2)*eta[1]*exp(z^2*(1-eta[2])/2)*(1-eta[2])*z*(eta[2]*(1-eta[1]) - (1-eta[1]))/((eta[2]^(1/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))^2)
	}
	deviance.mu <- function(z){
	               -2*log(((eta[2]^(1/2)*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))/(eta[2]^(1/2)*eta[1]*dnorm(0) + (1-eta[1])*dnorm(0))))
				   }
	tau <-  uniroot(function(x) v(x)*x^2 -1, lower=0, upper=35)$root
	deviance.phi <- function(z){
	               -2*log((abs(z/tau))*((eta[2]^(1/2)*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))/(eta[2]^(1/2)*eta[1]*dnorm(tau*sqrt(eta[2])) + (1-eta[1])*dnorm(tau))))
				   }
	cdfz <- function(z){
	        eta[1]*pnorm(z*sqrt(eta[2])) + (1-eta[1])*pnorm(z)
    }
	lpdf <- function(z,phi){
	       log(sqrt(eta[2])*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z)) -(1/2)*log(phi)

    }
	xix <- eta[1]/eta[2] + (1-eta[1])
 	dgd <- function(z){
	       (v(z)*z)^2*(sqrt(eta[2])*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))}
	dg <- 2*integrate(dgd,0,35)$value
 	fgd <- function(z){
	       (v(z))^2*z^4*(sqrt(eta[2])*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))}
	fg <- 2*integrate(fgd,0,35)$value
}

if(family=="Powerexp"){
	if(xi[1]<=-1  | xi[1]>=1) stop("the extra parameter must be within the interval (-1, 1)!!",call.=FALSE)
	family="Power Exponential"	
	kk <- xi[1]
	v <- function(z){
		 abs(z)^(-(2*kk/(1+kk)))/(1+kk)
	}
	vp <- function(z){
		  -2*kk*ifelse(z>=0,1,-1)*abs(z)^(-((3*kk + 1)/(1+kk)))/(1+kk)^2
	}
	dg <- 2^(1-kk)*gamma((3-kk)/2)/((1+kk)^2*gamma((1+kk)/2))
	fg <- (kk + 3)/(kk + 1)
	deviance.mu <- function(z){(abs(z))^(2/(1+kk))}
	deviance.phi <- function(z){(abs(z))^(2/(1+kk)) - (1+kk) - log(z^2/((1+kk)^(1+kk)))}
	cdfz <- function(z){
		pp <- 2/(kk+1)
	    sigmap <- (1+kk)^((kk+1)/2)
	    pnormp(z,mu=0,sigmap=sigmap,p=pp)
    }
	lpdf <- function(z,phi){
	       log(exp(-abs(z)^(2/(1+kk))/2)/(gamma(1 + (1+kk)/2)*2^(1 + (1+kk)/2))) -(1/2)*log(phi)

    }
	xix <- 2^(1+kk)*gamma(3*(1+kk)/2)/gamma((1+kk)/2)
}
if(family=="Sinh-normal"){
   	if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	alpha <- xi[1]
	v <- function(z){
		 4*sinh(z)*cosh(z)/(alpha^2*z) - tanh(z)/z
	}
	vp <- function(z){
		  ((cosh(z)*z-sinh(z))/z^2)*(4*cosh(z)/alpha^2 - 1/cosh(z)) + (sinh(z)/z)*(4*sinh(z)/alpha^2 + sinh(z)/(cosh(z)^2))
	}
	dg <- 2 + 4/(alpha^2) - (sqrt(2*pi)/alpha)*(1-2*(pnorm(sqrt(2)/alpha,mean=0,sd=sqrt(2)/2)-0.5))*exp(2/(alpha^2))
	dshn <- function(z){
	        2*cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)/(alpha*sqrt(2*pi))}
 	fgf <- function(z){
	       dshn(z)*(4*sinh(z)*cosh(z)/(alpha^2) - tanh(z))^2*z^2}
	fg <- 2*integrate(fgf,0,60)$value
    deviance.mu <- function(z){
				   if(alpha<=2) 4*(sinh(z))^2/alpha^2 - log((cosh(z))^2)
				   else{
				       z0 <- acosh(alpha/2)
					   2*log(dshn(z0)/dshn(z))
				   }
	}
	tau <- uniroot(function(x) 4*x*sinh(x)*cosh(x)/(alpha^2) - tanh(x)*x -1, lower=0, upper=50)$root
    deviance.phi <- function(z){
	            a <- 2*log(cosh(tau)*exp(-(2/alpha^2)*(sinh(tau))^2)) + log(tau^2)                                       
				b <- 2*log(cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)) + log(z^2)
				ifelse(a<b,0,a-b)
	}
	cdfz <- function(z){pnorm(2*sinh(z)/alpha)}
	lpdf <- function(z,phi){
	       log((2*cosh(sqrt(z^2))*exp(-2*sinh(sqrt(z^2))*sinh(sqrt(z^2))/alpha^2))) -(1/2)*log(2*pi*alpha^2) -(1/2)*log(phi)
    }
	fgf <- function(z){
	       dshn(z)*z^2}
	xix <- 2*integrate(fgf,0,20)$value
	fg2 <- fg
}	

if(family=="Sinh-t"){
   	if(xi[1]<=0 | xi[2]<=0) stop("the extra parameters must be positive!!",call.=FALSE)
	alpha <- xi[1]
	nu <- xi[2]	
	v <- function(z){
		 4*sinh(z)*cosh(z)*(nu+1)/(alpha^2*z*nu + 4*z*sinh(z)*sinh(z)) - tanh(z)/z
	}
	vp <- function(z){
		  ((cosh(z)*z - sinh(z))/z^2)*(4*cosh(z)*(nu+1)/(alpha^2*nu + 4*sinh(z)*sinh(z)) - 1/cosh(z)) +
		  (sinh(z)/z)*((alpha^2*nu + 4*sinh(z)*sinh(z))*4*(nu+1)*sinh(z) - 32*(nu+1)*cosh(z)*cosh(z)*sinh(z))/(alpha^2*nu + 4*sinh(z)*sinh(z))^2
	}
	dsht <- function(z){
	        cosh(z)*(1 + 4*(sinh(z))^2/(nu*alpha^2))^(-(nu+1)/2)}
 	dgd <- function(z){
	       dsht(z)*(4*sinh(z)*cosh(z)*(nu+1)/(alpha^2*nu + 4*sinh(z)*sinh(z)) - tanh(z))^2/(2*integrate(dsht,0,80)$value)}
	dg <- 2*integrate(dgd,0,60)$value
 	fgf <- function(z){
	       dsht(z)*(4*sinh(z)*cosh(z)*(nu+1)/(alpha^2*nu + 4*sinh(z)*sinh(z)) - tanh(z))^2*z^2/(2*integrate(dsht,0,80)$value)}
	fg <- 2*integrate(fgf,0,70)$value
    deviance.mu <- function(z){
				   if(alpha<=2*sqrt(1 + 1/nu)) (nu + 1)*log(1 + 4*(sinh(z))^2/(nu*alpha^2)) - log((cosh(z))^2)
				   else{
				       z0 <- acosh(sqrt(alpha^2/4 - 1/nu))
					   2*log(dsht(z0)/dsht(z))
				   }
				   }
	tau <- uniroot(function(x) 4*sinh(x)*cosh(x)*(nu+1)*x/(alpha^2*nu + 4*sinh(x)*sinh(x)) - tanh(x)*x -1, lower=0, upper=80)$root
    deviance.phi <- function(z){
	            a <- 2*log(cosh(tau)*(1 + 4*(sinh(tau))^2/(nu*alpha^2))^(-(nu+1)/2)) + log(tau^2)                                       
				b <- 2*log(cosh(z)*(1 + 4*(sinh(z))^2/(nu*alpha^2))^(-(nu+1)/2)) + log(z^2)                                       
				ifelse(a<b,0,a-b)
	}
	cdfz <- function(z){pt(2*sinh(z)/alpha,nu)}
	lpdf <- function(z,phi){
	       log(dsht(z)/(2*integrate(dsht,0,80)$value)) -(1/2)*log(phi)
    }
	fgf <- function(z){
	       dsht(z)*z^2/(2*integrate(dsht,0,80)$value)}
	xix <- 2*integrate(fgf,0,50)$value
	dshn <- function(z){
	2*cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)/(alpha*sqrt(2*pi))}
 	fgf <- function(z){
	dshn(z)*(4*sinh(z)*cosh(z)/(alpha^2) - tanh(z))^2*z^2}
	fg2 <- 2*integrate(fgf,0,60)$value
}

if(family=="Hyperbolic"){
    if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	family="Symmetric Hyperbolic"	
	nu <- xi[1]
	v <- function(z){
		 nu/sqrt(1 + z^2)
	}
	vp <- function(z){
		  -nu*z/((1 + z^2)^(3/2))
	}
	dh <- function(z){
	      exp(-nu*sqrt(1+z^2))}
 	dgd <- function(z){
	       dh(z)*(nu*z/sqrt(1 + z^2))^2/(2*integrate(dh,0,Inf)$value)}
	dg <- 2*integrate(dgd,0,Inf)$value
 	fgf <- function(z){
	       dh(z)*(nu*z^2/sqrt(1 + z^2))^2/(2*integrate(dh,0,Inf)$value)}
	fg <- 2*integrate(fgf,0,Inf)$value
    deviance.mu <- function(z){2*nu*(sqrt(1+z^2)-1)}
	tau <- sqrt((1 + sqrt(1 + 4*nu^2))/(2*nu^2))
    deviance.phi <- function(z){2*nu*(sqrt(1+z^2) - sqrt(1 + tau^2)) - log(z^2/tau^2)}	
	cdfz <- function(z){temporal <- matrix(0,length(z),1)
	                    for(gg in 1:length(z)){
                       	    temporal[gg] <- integrate(dh,-Inf,z[gg])$value
						}
						temporal/(2*integrate(dh,0,Inf)$value)
    }
	lpdf <- function(z,phi){
	       log(dh(z)/(2*integrate(dh,0,Inf)$value)) -(1/2)*log(phi)
    }
	fgf <- function(z){
	       dh(z)*z^2/(2*integrate(dh,0,Inf)$value)}
	xix <- 2*integrate(fgf,0,Inf)$value
}

if(family=="Slash"){
   	if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	nu <- xi[1]
	G <- function(a,x){
	     gamma(a)*gamma_inc_P(a,x)/(x^a)
	}	 
	v <- function(z){
		 G(nu+3/2,z^2/2)/G(nu+1/2,z^2/2)
	}
	ds <- function(z){G(nu+1/2,z^2/2)}
	gdg <- function(z){ds(z)*(v(z))^2*z^2/(2*integrate(ds,0,Inf)$value)}
	dg <- 2*integrate(gdg,0,Inf)$value
	gfg <- function(z){ds(z)*(v(z))^2*z^4/(2*integrate(ds,0,Inf)$value)}
	fg <- 2*integrate(gfg,0,Inf)$value
    deviance.mu <- function(z){2*log(2/(2*nu+1))-2*log(G(nu+1/2,z^2/2))}
	tau <- uniroot(function(x) v(x)*x^2 -1, lower=0.0001, upper=1000)$root
    deviance.phi <- function(z){
	            a <- 2*log(G(nu+1/2,tau^2/2)) + log(tau^2)                                       
				b <- 2*log(G(nu+1/2,z^2/2)) + log(z^2)                                                                              
				ifelse(a<b,0,a-b)
	}
	cdfz <- function(z){temporal <- matrix(0,length(z),1)
	                    for(gg in 1:length(z)){
                       	    temporal[gg] <- integrate(ds,-Inf,z[gg])$value
						}
						temporal/(2*integrate(ds,0,Inf)$value)
    }
	lpdf <- function(z,phi){
	       log(ds(z)/(2*integrate(ds,0,Inf)$value)) -(1/2)*log(phi)
    }
	gfg <- function(z){ds(z)*z^2/(2*integrate(ds,0,Inf)$value)}
	xix <- 2*integrate(gfg,0,60)$value
}


if(p>0){
	if(l>0){
	objeto <- list(epsilon=epsilon,maxiter=maxiter,p=p,l=l,q=q,qm=qm,mu=mu,GradD=GradD,y=y,W=W,v=v,dg=dg,fg=fg,
	              subset=subset,family=family)
	}
	else{
	objeto <- list(epsilon=epsilon,maxiter=maxiter,p=p,l=l,q=q,qm=qm,mu=mu,GradD=GradD,y=y,W=matrix(1,n,1),v=v,dg=dg,fg=fg,
			      subset=subset,family=family)
	}
}else{if(l>0){
	objeto <- list(epsilon=epsilon,maxiter=maxiter,p=p,l=l,q=q,qm=qm,y=y,W=W,v=v,dg=dg,fg=fg,
	              subset=subset,family=family)
	}
	else{
	objeto <- list(epsilon=epsilon,maxiter=maxiter,p=p,l=l,q=q,qm=qm,y=y,W=matrix(1,n,1),v=v,dg=dg,fg=fg,
			      subset=subset,family=family)
	}
}


if(qm>0){
if(missingArg(lambda.mu)){
  if(p>0) {
    X_au <- cbind(X,1,xxm,xxm^2)
	b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%y
  	yres <- y - X%*%(b_au[1:p])
  	yres2 <- (y - X_au%*%b_au)^2
  }
  if(p==0) {
    X_au <- cbind(1,xxm,xxm^2)
	b_au <- solve(t(X_au)%*%X_au)%*%t(X_au)%*%y
  	yres <- y
  	yres2 <- (y - X_au%*%b_au)^2
  }
  lambda.mu_est <- lambda.hat(yres,xxm,start.lambda.mu)
  lambda.mu <- lambda.mu_est$lambda_hat/median(yres2)
  l.mu <- lambda.mu*median(yres2)  
}else{l.mu <- lambda.mu}

  ssm <- splinek(xxm)
  Nm <- ssm$N
  Mm <- ssm$K
  objeto$Nm <- Nm
  objeto$Mm <- Mm
  objeto$lambda.mu <- lambda.mu
  Hm <- Nm%*%solve(t(Nm)%*%Nm + l.mu*Mm)%*%t(Nm)
}


if(q>0){
  if(missingArg(lambda.phi)){
	if(qm>0){
	  objeto$q <- 0
	  objeto$W <- matrix(1,n,1)
	  objeto$l <- 1
	  objeto$K_psi2 <- 2/n
	  objeto$W_bar <- objeto$W
	  if(p>0){
		  if(family=="Sinh-t" | family=="Sinh-normal"){
		       objeto$xi <- xi
			   objeto$K_psi2 <- 4/((fg2-1)*n)
		       vP <- itpE2(c(beta0, rep(0,qm),log(mean((y-mu(beta0))^2))),objeto)
		  }else{
		       vP <- itpE(c(beta0, rep(0,qm),log(mean((y-mu(beta0))^2))),objeto)
		  }
		  phi.sint <- log((y - mu(vP$theta[1:p])-Nm%*%vP$theta[(p+1):(p+qm)])^2)
	  }else{
		  if(family=="Sinh-t" | family=="Sinh-normal"){
		       objeto$xi <- xi
			   objeto$K_psi2 <- 4/((fg2-1)*n)
		       vP <- itpE2(c(rep(mean(y),qm),log(mean((y)^2))),objeto)
		  }else{
		       vP <- itpE(c(rep(mean(y),qm),log(mean((y)^2))),objeto)
		  }
		  phi.sint <- log((y -Nm%*%vP$theta[(p+1):(p+qm)])^2)
	   }
	  objeto$q <- q
	  if(l>0) objeto$W <- W
	  objeto$l <- l  
	  if(l>0){
	    W_au <- cbind(W,1,xx,xx^2)
	    phi.sint <- phi.sint - W%*%(solve(t(W_au)%*%W_au)%*%t(W_au)%*%phi.sint)[1:l]
	  }
	  lambda.phi_est <- lambda.hat(phi.sint,xx,start.lambda.phi)
	  l.phi <- lambda.phi_est$lambda_hat
	}else{
  	  objeto$q <- 0
	  objeto$W <- matrix(1,n,1)
	  objeto$l <- 1
	  objeto$K_psi2 <- 2/n
	  objeto$W_bar <- objeto$W
		  if(family=="Sinh-t" | family=="Sinh-normal"){
		       objeto$xi <- xi
			   objeto$K_psi2 <- 4/((fg2-1)*n)
		       vP <- itpE2(c(beta0,log(mean((y-mu(beta0))^2))),objeto)
		  }else{
		       vP <- itpE(c(beta0,log(mean((y-mu(beta0))^2))),objeto)
		  }
	  objeto$q <- q
	  if(l>0) objeto$W <- W
	  objeto$l <- l  
	  phi.sint <- log((y - mu(vP$theta[1:p]))^2)
	  if(l>0){
	    W_au <- cbind(W,1,xx,xx^2)
	    phi.sint <- phi.sint - W%*%(solve(t(W_au)%*%W_au)%*%t(W_au)%*%phi.sint)[1:l]
	  }
	  lambda.phi_est <- lambda.hat(phi.sint,xx,start.lambda.phi)
	  l.phi <- lambda.phi_est$lambda_hat/2
	}
  }else{l.phi <- lambda.phi}
  ss <- splinek(xx)
  N <- ss$N
  M <- ss$K
  objeto$N <- N
  objeto$M <- M
  objeto$lambda.phi <- l.phi
  gle <- sum(diag(N%*%solve(t(N)%*%N + l.phi*M)%*%t(N)))
}

if(p>0){
	if(qm==0){
		if(q>0){
		  if(missingArg(lambda.phi)){
			if(l==0) {theta0 <- c(vP$theta[1:p], rep(vP$theta[(p+1)],q))}
			if(l>0){theta0 <- c(vP$theta[1:p], rep(0,l), rep(vP$theta[(p+1)],q))}
		  }
		  else{
			if(l==0) {theta0 <- c(beta0, rep(mean(log((y-mu(beta0))^2)),q))}
			if(l>0) {theta0 <- c(beta0, rep(0,l), rep(mean(log((y-mu(beta0))^2)),q))}
		  }
		}													
		if(q==0){
		  if(missingArg(formula.phi)){
		        theta0 <- c(beta0, mean(log((y-mu(beta0))^2/ifelse(is.null(xix),1,xix))))
		  }
		  else{
		     aa <- solve(t(W)%*%W)%*%t(W)%*%(log((y-mu(beta0))^2/ifelse(is.null(xix),1,xix)))
			 theta0 <- c(beta0, aa)
		  }
		}	   
	}
	if(qm>0){
		if(q>0){
		  if(missingArg(lambda.phi)){
			if(l==0) {theta0 <- c(vP$theta[1:(p+qm)], rep(vP$theta[(p+qm+1)],q))}
			if(l>0){theta0 <- c(vP$theta[1:(p+qm)], rep(0,l), rep(vP$theta[(p+qm+1)],q))}
		  }
		  else{
			if(l==0) {theta0 <- c(beta0, rep(0,qm),rep(mean(log((y-mu(beta0))^2/ifelse(is.null(xix),1,xix))),q))}
			else{theta0 <- c(beta0, rep(0,qm),rep(0,l),rep(mean(log((y-mu(beta0))^2/ifelse(is.null(xix),1,xix))),q))}
		  }
		}
		if(q==0){
		  yres <- y - X%*%solve(t(X)%*%X)%*%t(X)%*%y
		  if(missingArg(formula.phi)){
		        theta0 <- c(beta0, rep(0,qm),mean(log((y-mu(beta0)-Hm%*%yres)^2/ifelse(is.null(xix),1,xix))))
		  }
		  else{
		     aa <- solve(t(W)%*%W)%*%t(W)%*%(log((y-mu(beta0)-Hm%*%yres)^2/ifelse(is.null(xix),1,xix)))
			 theta0 <- c(beta0, rep(0,qm), aa)
		  }
		}	   
	}
}else{
	if(qm>0){
		if(q>0){
		  if(missingArg(lambda.phi)){
			if(l==0) {theta0 <- c(vP$theta[1:qm], rep(vP$theta[(qm+1)],q))}
			else{theta0 <- c(vP$theta[1:qm], rep(0,l), rep(vP$theta[(qm+1)],q))}
		  }
		  else{
			if(l==0) {theta0 <- c(solve(t(Nm)%*%Nm + l.mu*Mm)%*%t(Nm)%*%y, rep(mean(log((y-Hm%*%y)^2/ifelse(is.null(xix),1,xix))),q))}
			else{theta0 <- c(solve(t(Nm)%*%Nm + l.mu*Mm)%*%t(Nm)%*%y, rep(0,l), rep(mean(log((y-Hm%*%y)^2/ifelse(is.null(xix),1,xix))),q))}
		  }
		}
		if(q==0){
		  if(missingArg(formula.phi)){
		        theta0 <- c(solve(t(Nm)%*%Nm + l.mu*Mm)%*%t(Nm)%*%y,mean(log((y-Hm%*%y)^2/ifelse(is.null(xix),1,xix))))
		  }
		  else{
		     aa <- solve(t(W)%*%W)%*%t(W)%*%log((y-Hm%*%y)^2/ifelse(is.null(xix),1,xix))
			 if(l>1) theta0 <- c(solve(t(Nm)%*%Nm + l.mu*Mm)%*%t(Nm)%*%y, aa)
		  }
		}	   
	}
}

if(q>0) lambda.phi <- l.phi


	if(l>0){
	  if(q>0){
			W_bar <- cbind(W,objeto$N)
			K_psi <- t(W_bar)%*%W_bar
			if(family=="Sinh-t" | family=="Sinh-normal"){
				K_psi[(l+1):(l+q),(l+1):(l+q)] <- K_psi[(l+1):(l+q),(l+1):(l+q)] + (4*objeto$lambda.phi/(fg2-1))*objeto$M
			    K_psi2 <- solve(((fg2-1)/4)*K_psi)
            }else{
				K_psi[(l+1):(l+q),(l+1):(l+q)] <- K_psi[(l+1):(l+q),(l+1):(l+q)] + (2*objeto$lambda.phi)*objeto$M
			    K_psi2 <- solve((1/2)*K_psi)
			}   
	  }
	  if(q==0){
	    W_bar <- W
		K_psi <- t(W_bar)%*%W_bar	
	    if(family=="Sinh-t" | family=="Sinh-normal"){
		        K_psi2 <- solve(((fg2-1)/4)*K_psi)
		}else{
		     K_psi2 <- solve((1/2)*K_psi)		     
		}
	  }
	}
	if(l==0){
	  if(q>0){
			W_bar <- objeto$N
			if(family=="Sinh-t" | family=="Sinh-normal"){
				K_psi <- t(W_bar)%*%W_bar + (4*objeto$lambda.phi/(fg2-1))*objeto$M
			    K_psi2 <- solve(((fg2-1)/4)*K_psi)
            }else{
				K_psi <- t(W_bar)%*%W_bar + (2*objeto$lambda.phi)*objeto$M
			    K_psi2 <- solve((1/2)*K_psi)
			}   
	  }
    }

objeto$K_psi2 <- K_psi2
objeto$W_bar <- W_bar
if(family=="Sinh-t" | family=="Sinh-normal"){
       objeto$xi <- xi
	  ajuste <- itpE2(theta0,objeto)
}else{
      ajuste <- itpE(theta0,objeto)
}   



theta <- ajuste$theta 

if(p>0){
  mues <- mu(theta[1:p])[subset,]
  if(qm==0){
	if(l>0){
	  if(q>0){
	    	gammav <- theta[(p+1):(p+l)]
		    h <- theta[(p+l+1):length(theta)]
			phi.es <- exp(W%*%gammav + objeto$N%*%h)
			z <- (y-mues)/sqrt(phi.es)			
			psi <- theta[(p+1):length(theta)]
			W_bar <- cbind(W,objeto$N)
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
			K_psi <- t(W_bar)%*%W_bar	
			K_psi[(l+1):length(psi),(l+1):length(psi)] <- K_psi[(l+1):length(psi),(l+1):length(psi)] + (4*objeto$lambda.phi/(fg-1))*objeto$M
	    	K_psi <- ((fg-1)/4)*K_psi
	    	K_psi2 <- solve(K_psi)
		
	  }
	  if(q==0){
    	     gammav <- theta[(p+1):(p+l)]
   			 psi <- theta[(p+1):length(theta)]
		     phi.es <- exp(W%*%gammav)
			 z <- (y-mues)/sqrt(phi.es)			 
			 W_bar <- W
			 U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			 K_psi <- t(W_bar)%*%W_bar	
	    	 K_psi <- ((fg-1)/4)*K_psi
	    	 K_psi2 <- solve(K_psi)
	  }
	}
	if(l==0){
	  if(q>0){
			h <- theta[(p+1):length(theta)]
			psi <- theta[(p+1):length(theta)]
		    phi.es <- exp(objeto$N%*%h)
			z <- (y-mues)/sqrt(phi.es)
			W_bar <- objeto$N
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			K_psi <- t(W_bar)%*%W_bar + (4*objeto$lambda.phi/(fg-1))*objeto$M
		    K_psi <- ((fg-1)/4)*K_psi
		    K_psi2 <- solve(K_psi)
	  }
    }
    Kinver <- K_psi2
    D <- GradD(theta[1:p])[objeto$subset,]
	D_bar <- D
	D2 <- D_bar*matrix(1/phi.es,length(phi.es),p+qm)
	se.beta <- sqrt(diag(solve(dg*t(D_bar)%*%D2)))
	se.psi <- sqrt(diag(((fg-1)/4)*Kinver%*%(t(W_bar)%*%W_bar)%*%Kinver))
	escore <- t(D_bar)%*%(v(z)*((y-mues)/phi.es))
	escore <- rbind(escore,U_psi)

  }

  if(qm>0){
	if(l>0){
	  if(q>0){
	  		hm <- theta[(p+1):(p+qm)]
	    	gammav <- theta[(p+qm+1):(p+qm+l)]
		    h <- theta[(p+qm+l+1):length(theta)]
			phi.es <- exp(W%*%gammav + objeto$N%*%h)
			z <- (y-mues-Nm%*%hm)/sqrt(phi.es)			
			psi <- theta[(p+qm+1):length(theta)]
			W_bar <- cbind(W,objeto$N)
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
			K_psi <- t(W_bar)%*%W_bar	
			K_psi[(l+1):length(psi),(l+1):length(psi)] <- K_psi[(l+1):length(psi),(l+1):length(psi)] + (4*objeto$lambda.phi/(fg-1))*objeto$M
	    	K_psi <- ((fg-1)/4)*K_psi
	    	K_psi2 <- solve(K_psi)
	  }
	  if(q==0){
	  		 hm <- theta[(p+1):(p+qm)]
    	     gammav <- theta[(p+qm+1):(p+qm+l)]
   			 psi <- theta[(p+qm+1):length(theta)]
		     phi.es <- exp(W%*%gammav)
			 z <- (y-mues-Nm%*%hm)/sqrt(phi.es)			 
			 W_bar <- W
			 U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			 K_psi <- t(W_bar)%*%W_bar	
	    	 K_psi <- ((fg-1)/4)*K_psi
	    	 K_psi2 <- solve(K_psi)
	  }
	}
	if(l==0){
	  if(q>0){
	  		hm <- theta[(p+1):(p+qm)]
			h <- theta[(p+qm+1):length(theta)]
			psi <- theta[(p+qm+1):length(theta)]
		    phi.es <- exp(objeto$N%*%h)
			z <- (y-mues-Nm%*%hm)/sqrt(phi.es)
			W_bar <- objeto$N
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			K_psi <- t(W_bar)%*%W_bar + (4*objeto$lambda.phi/(fg-1))*objeto$M
		    K_psi <- ((fg-1)/4)*K_psi
		    K_psi2 <- solve(K_psi)
	  }
    }
	Kinver <- K_psi2
	D <- GradD(theta[1:p])[objeto$subset,]
	D_bar <- cbind(D,Nm)
	Kpqm <- matrix(0,p+qm,p+qm)
	Kpqm[(p+1):(p+qm),(p+1):(p+qm)] <- lambda.mu*Mm
	D2 <- D_bar*matrix(1/phi.es,length(phi.es),p+qm)
	se.beta <- sqrt(diag(solve(dg*t(D_bar)%*%D2 + Kpqm)%*%(dg*t(D_bar)%*%D2)%*%solve(dg*t(D_bar)%*%D2 + Kpqm)))
	se.psi <- sqrt(diag(((fg-1)/4)*Kinver%*%(t(W_bar)%*%W_bar)%*%Kinver))
	escore <- t(D_bar)%*%(v(z)*((y-mues-Nm%*%hm)/phi.es))  - rbind(matrix(0,p,1),lambda.mu*Mm%*%hm)
	escore <- rbind(escore,U_psi)
	Nme <- Nm*matrix(1/sqrt(phi.es),length(phi.es),ncol(Nm))
	Hm <- Nm%*%solve(t(Nme)%*%Nme + lambda.mu*Mm)%*%t(Nme*matrix(1/sqrt(phi.es),length(phi.es),ncol(Nm)))
    gle.mu <- sum(diag(Hm))

  }
}
if(p==0){
  if(qm>0){
	if(l>0){
	  if(q>0){
	  		hm <- theta[(p+1):(p+qm)]
	    	gammav <- theta[(p+qm+1):(p+qm+l)]
		    h <- theta[(p+qm+l+1):length(theta)]
			phi.es <- exp(W%*%gammav + objeto$N%*%h)
			z <- (y-Nm%*%hm)/sqrt(phi.es)			
			psi <- theta[(p+qm+1):length(theta)]
			W_bar <- cbind(W,objeto$N)
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
			K_psi <- t(W_bar)%*%W_bar	
			K_psi[(l+1):length(psi),(l+1):length(psi)] <- K_psi[(l+1):length(psi),(l+1):length(psi)] + (4*objeto$lambda.phi/(fg-1))*objeto$M
	    	K_psi <- ((fg-1)/4)*K_psi
	    	K_psi2 <- solve(K_psi)
	  }
	  if(q==0){
	  		 hm <- theta[(p+1):(p+qm)]
    	     gammav <- theta[(p+qm+1):(p+qm+l)]
   			 psi <- theta[(p+qm+1):length(theta)]
		     phi.es <- exp(W%*%gammav)
			 z <- (y-Nm%*%hm)/sqrt(phi.es)			 
			 W_bar <- W
			 U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			 K_psi <- t(W_bar)%*%W_bar	
	    	 K_psi <- ((fg-1)/4)*K_psi
	    	 K_psi2 <- solve(K_psi)
	  }
	}
	if(l==0){
	  if(q>0){
	  		hm <- theta[(p+1):(p+qm)]
			h <- theta[(p+qm+1):length(theta)]
			psi <- theta[(p+qm+1):length(theta)]
		    phi.es <- exp(objeto$N%*%h)
			z <- (y-Nm%*%hm)/sqrt(phi.es)
			W_bar <- objeto$N
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			K_psi <- t(W_bar)%*%W_bar + (4*objeto$lambda.phi/(fg-1))*objeto$M
		    K_psi <- ((fg-1)/4)*K_psi
		    K_psi2 <- solve(K_psi)
			
	  }
    }
	Kinver <- K_psi2
	D_bar <- Nm
	D2 <- D_bar*matrix(1/phi.es,length(phi.es),p+qm)
	se.beta <- sqrt(diag(solve(dg*t(D_bar)%*%D2 + lambda.mu*Mm)%*%(dg*t(D_bar)%*%D2)%*%solve(dg*t(D_bar)%*%D2 + lambda.mu*Mm)))
	se.psi <- sqrt(diag(((fg-1)/4)*Kinver%*%(t(W_bar)%*%W_bar)%*%Kinver))
	escore <- t(D_bar)%*%(v(z)*((y-Nm%*%hm)/phi.es))  - lambda.mu*Mm%*%hm
	escore <- rbind(escore,U_psi)
	Nme <- Nm*matrix(1/sqrt(phi.es),length(phi.es),ncol(Nm))
	Hm <- Nm%*%solve(t(Nme)%*%Nme + lambda.mu*Mm)%*%t(Nme*matrix(1/sqrt(phi.es),length(phi.es),ncol(Nm)))
    gle.mu <- sum(diag(Hm))

  }
}

if(!missingArg(local.influence)){
    if(family=="Slash"){
	   vp <- grad(v,z)
	   Da <- (vp*z + v(z))
	   Dh <- (vp*z^2 + 2*z*v(z))/2
	   Dc <- (vp*z^3 + 2*z^2*v(z))/4
	}else{
	Da <- (vp(z)*z + v(z))
	Dh <- (vp(z)*z^2 + 2*z*v(z))/2
	Dc <- (vp(z)*z^3 + 2*z^2*v(z))/4
	}
	Lbb <- -t(D_bar)%*%(matrix(Da/phi.es,length(phi.es),ncol(D_bar))*D_bar)
	Lgg <- -t(W_bar)%*%(matrix(Dc,nrow(W_bar),ncol(W_bar))*W_bar)
	Lbg <- -t(D_bar)%*%(matrix(Dh/sqrt(phi.es),nrow(W_bar),ncol(W_bar))*W_bar)
	if(qm>0){
	  Lbb[(p+1):(p+qm),(p+1):(p+qm)] <-   Lbb[(p+1):(p+qm),(p+1):(p+qm)] - lambda.mu*Mm
	}
	if(q>0){
	  Lgg[(l+1):ncol(W_bar),(l+1):ncol(W_bar)] <-   Lgg[(l+1):ncol(W_bar),(l+1):ncol(W_bar)] - lambda.phi*M
	}
	Ltt <- matrix(0,p+qm+l+q,p+qm+l+q)
	Ltt2 <- matrix(0,p+qm+l+q,p+qm+l+q)
	Ltt[1:(p+qm),1:(p+qm)] <- Lbb
	Ltt[1:(p+qm),(p+qm+1):(p+qm+l+q)]  <- Lbg
	Ltt[(p+qm+1):(p+qm+l+q),1:(p+qm)]  <- t(Lbg)
	Ltt[(p+qm+1):(p+qm+l+q),(p+qm+1):(p+qm+l+q)] <- Lgg
	Ltt2[(p+qm+1):(p+qm+l+q),(p+qm+1):(p+qm+l+q)] <- solve(Lgg)

	delta1 <- t(D_bar*matrix(v(z)*z/sqrt(phi.es),length(phi.es),p+qm))
	delta2 <- (1/2)*t(matrix(v(z)*z^2-1,nrow(W_bar),ncol(W_bar))*W_bar)
	delta <- rbind(delta1,delta2)
	
	tl <- t(delta)%*%(solve(Ltt)-Ltt2)%*%delta
	tl <- tl/sqrt(sum(diag(t(tl)%*%tl)))
	cw <- abs(eigen(tl,symmetric=TRUE)$vector[,length(y)])
	cw2 <- abs(diag(tl))

	delta1 <- t(D_bar*matrix(Da/phi.es,length(phi.es),p+qm))
	delta2 <- t(matrix(Dh/sqrt(phi.es),nrow(W_bar),ncol(W_bar))*W_bar)
	delta <- rbind(delta1,delta2)
	
	tl2 <- t(delta)%*%(solve(Ltt)-Ltt2)%*%delta
	tl2 <- tl2/sqrt(sum(diag(t(tl2)%*%tl2)))
	pr <- abs(eigen(tl2,symmetric=T)$vector[,length(y)])
	pr2 <- abs(diag(tl2))
	
	if(p>0){
	 if(qm>0){
	   if(l>0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
	   }
	   if(l==0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
	   }
	 }
	 if(qm==0){
	   if(l>0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
	   }
	   if(l==0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
	   }
	 }
	}
	if(p==0){
	 if(qm>0){
	   if(l>0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas2=colnames(W),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas2=colnames(W),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
	   }
	   if(l==0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),xix=xix,cw=cbind(cw,cw2),pr=cbind(pr,pr2),call="")
		 }
	   }
	 }
	}
}
else{
	if(p>0){
	 if(qm>0){
	   if(l>0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,call="")
		 }
	   }
	   if(l==0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,call="")
		 }
	   }
	 }
	 if(qm==0){
	   if(l>0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),filas2=colnames(W),xix=xix,call="")
		 }
	   }
	   if(l==0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas=colnames(X),xix=xix,call="")
		 }
	   }
	 }
	}
	if(p==0){
	 if(qm>0){
	   if(l>0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas2=colnames(W),xix=xix,call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),filas2=colnames(W),xix=xix,call="")
		 }
	   }
	   if(l==0){
	     if(q>0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.phi=gle,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,lambda.phi=lambda.phi,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),xix=xix,call="")
		 }
		 if(q==0){
			 salida <- list(p=p,l=l,q=q,qm=qm,gle.mu=gle.mu,coefs.mu=ajuste$theta[1:(p+qm)],se.mu=se.beta,coefs.phi=ajuste$theta[(p+qm+1):(p+qm+l+q)],se.phi=se.psi,
			 z.hat=z,xi=xi,lambda.mu=lambda.mu,family=family,v=v(z),dg=dg,fg=fg,cdfz=cdfz(z),mu.fitted=y-z*sqrt(phi.es), phi.fitted=phi.es,deviance.mu=deviance.mu(z),
			 deviance.phi=deviance.phi(z),lpdf=lpdf(z,phi.es),xix=xix,call="")
		 }
	   }
	 }
	}
}

 class(salida) <- "ssym"
 salida$call <- match.call()

salida		
}
