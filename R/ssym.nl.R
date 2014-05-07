ssym.nl <-
function(response, mu, start, formula.phi, ncs, start.lambda, lambda, family, xi, epsilon, maxiter, subset, local.influence){

if(family!="Normal" & family!="Slash" & family!="Hyperbolic" & family!="Sinh-t" & family!="Sinh-normal" & family!="Contnormal" & family!="Powerexp" & family!="Student")
stop("family of distributions specified by the user is not supported!!",call.=FALSE)

if(family=="Slash" | family=="Hyperbolic" | family=="Sinh-t" | family=="Sinh-normal" | family=="Contnormal" | family=="Powerexp" | family=="Student"){
  if(missingArg(xi)) stop("for the family of distributions specified by the user an extra parameter is required!!", call.=FALSE) 
}

beta0 <- start
p <- length(start)

GradD <- function(vP){jacobian(mu,vP[1:p])}
if(missingArg(epsilon)) epsilon <- 0.0000001
if(missingArg(maxiter)) maxiter <- 1000
if(missingArg(start.lambda)) start.lambda <- 1
if(missingArg(subset)) subset <- 1:length(response)

y <- response[subset]
n <- length(y)


if(missingArg(ncs)){
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
   xx <- ncs[subset]
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


if(l>0){
objeto <- list(epsilon=epsilon,maxiter=maxiter,p=p,l=l,q=q,qm=0,mu=mu,GradD=GradD,y=y,W=W,v=v,dg=dg,fg=fg,N=matrix(0,3,1),
M=matrix(0,3,1),lambda=1,subset=subset,family=family)
}
else{
objeto <- list(epsilon=epsilon,maxiter=maxiter,p=p,l=l,q=q,qm=0,mu=mu,GradD=GradD,y=y,W=matrix(1,n,1),v=v,dg=dg,fg=fg,
N=matrix(0,3,1),M=matrix(0,3,1),lambda=1,subset=subset,family=family)
}

if(q>0){
if(missingArg(lambda)){
  objeto$q <- 0
  objeto$W <- matrix(1,n,1)
  objeto$l <- 1
  objeto$K_psi2 <- 2/n
  objeto$W_bar <- objeto$W
  if(family=="Sinh-t" | family=="Sinh-normal"){
       objeto$xi <- xi
	   objeto$K_psi2 <- 4/((fg2-1)*n)
       vP <- itpE2(c(beta0, log(mean((y-mu(beta0)[subset])^2))),objeto)
  }else{
       vP <- itpE(c(beta0, log(mean((y-mu(beta0)[subset])^2))),objeto)
  }
  objeto$q <- q
  if(l>0) objeto$W <- W
  objeto$l <- l  
  phi.sint <- log((y - mu(vP$theta)[subset])^2)
  if(l>0){
    W_au <- cbind(W,1,xx,xx^2)
    phi.sint <- phi.sint - W%*%(solve(t(W_au)%*%W_au)%*%t(W_au)%*%phi.sint)[1:l]
  }
  lambda_est <- lambda.hat(phi.sint,xx,start.lambda)
  lambd <- lambda_est$lambda_hat/2
}else{lambd <- lambda}

  ss <- splinek(xx)
  N <- ss$N
  M <- ss$K
  objeto$N <- N
  objeto$M <- M
  objeto$lambda.phi <- lambd
  gle <- sum(diag(N%*%solve(t(N)%*%N + lambd*M)%*%t(N)))
}
		if(q>0){
		  if(missingArg(lambda)){
			if(l==0) {theta0 <- c(vP$theta[1:p], rep(vP$theta[(p+1)],q))}
			if(l>0){theta0 <- c(vP$theta[1:p], rep(0,l), rep(vP$theta[(p+1)],q))}
		  }
		  else{
			if(l==0) {theta0 <- c(beta0, rep(mean(log((y-mu(beta0)[subset])^2/ifelse(is.null(xix),1,xix))),q))}
			if(l>0) {theta0 <- c(beta0, rep(0,l), rep(mean(log((y-mu(beta0)[subset])^2/ifelse(is.null(xix),1,xix))),q))}
		  }
		}
		if(q==0){
		  if(missingArg(formula.phi)){
		        theta0 <- c(beta0, mean(log((y-mu(beta0)[subset])^2/ifelse(is.null(xix),1,xix))))
		  }
		  else{
		     aa <- solve(t(W)%*%W)%*%t(W)%*%(log((y-mu(beta0)[subset])^2/ifelse(is.null(xix),1,xix)))
			 theta0 <- c(beta0, aa)
		  }
		}	   


	  if(q>0){ lambda <- lambd}

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
beta <- ajuste$theta[1:p]
mues <- mu(beta)[objeto$subset]

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
iter <- ajuste$iter

z.hat <- (y-mues)/sqrt(phi.es)
D <- GradD(beta)[objeto$subset,]
D2 <- D*matrix(1/phi.es,length(phi.es),p)
vcov.beta <- solve(dg*t(D)%*%D2)
se.psi <- sqrt(diag(((fg-1)/4)*Kinver%*%(t(W_bar)%*%W_bar)%*%Kinver))
escore <- t(D)%*%(v(z.hat)*((y-mues)/phi.es))
escore <- rbind(escore,U_psi)

D3 <- D*matrix(1/sqrt(phi.es),length(phi.es),p)
h <- dg*diag(D3%*%vcov.beta%*%t(D3))

if(!missingArg(local.influence)){
   Hessian <- array(0,c(p,p,n))
   for(ij in 1:n){
   		mu_ij <- function(vP){mu(vP)[ij]}
   		Hessian[,,ij] <- hessian(func=mu_ij, x=beta)
   }
	Dbb <- Hessian[,,objeto$subset]
    if(family=="Slash"){
	   vp <- grad(v,z.hat)
	   Da <- (vp*z.hat + v(z.hat))
	   Dh <- (vp*z.hat^2 + 2*z.hat*v(z.hat))/2
	   Dc <- (vp*z.hat^3 + 2*z.hat^2*v(z.hat))/4
	}else{
	Da <- (vp(z.hat)*z.hat + v(z.hat))
	Dh <- (vp(z.hat)*z.hat^2 + 2*z.hat*v(z.hat))/2
	Dc <- (vp(z.hat)*z.hat^3 + 2*z.hat^2*v(z.hat))/4
	}

	splb <- matrix(0,p,p)
	for(i in 1:length(z.hat)){
		splb <- splb + (z.hat[i]*v(z.hat)[i]/sqrt(phi.es[i]))*Dbb[,,i]
	}
	
	Lbb <- -t(D)%*%(matrix(Da/phi.es,length(phi.es),p)*D) + splb
	Lgg <- -t(W_bar)%*%(matrix(Dc,nrow(W_bar),ncol(W_bar))*W_bar)
	Lbg <- -t(D)%*%(matrix(Dh/sqrt(phi.es),nrow(W_bar),ncol(W_bar))*W_bar)
	if(q>0){
	  Lgg[(l+1):ncol(W_bar),(l+1):ncol(W_bar)] <-   Lgg[(l+1):ncol(W_bar),(l+1):ncol(W_bar)] - objeto$lambda*objeto$M
	}
	Ltt <- matrix(0,p+ncol(W_bar),p+ncol(W_bar))
	Ltt2 <- matrix(0,p+ncol(W_bar),p+ncol(W_bar))
	Ltt[1:p,1:p] <- Lbb
	Ltt[1:p,(p+1):ncol(Ltt)]  <- Lbg
	Ltt[(p+1):ncol(Ltt),1:p]  <- t(Lbg)
	Ltt[(p+1):ncol(Ltt),(p+1):ncol(Ltt)] <- Lgg
	Ltt2[(p+1):ncol(Ltt),(p+1):ncol(Ltt)] <- solve(Lgg)
	
	delta1 <- t(D*matrix(v(z.hat)*z.hat/sqrt(phi.es),length(phi.es),p))
	delta2 <- (1/2)*t(matrix(v(z.hat)*z.hat^2-1,nrow(W_bar),ncol(W_bar))*W_bar)
	delta <- rbind(delta1,delta2)
	
	tl <- t(delta)%*%(solve(Ltt)-Ltt2)%*%delta
	tl <- tl/sqrt(sum(diag(t(tl)%*%tl)))
	cw <- abs(eigen(tl,symmetric=TRUE)$vector[,length(y)])
	cw2 <- abs(diag(tl))

	delta1 <- t(D*matrix(Da/phi.es,length(phi.es),p))
	delta2 <- t(matrix(Dh/sqrt(phi.es),nrow(W_bar),ncol(W_bar))*W_bar)
	delta <- rbind(delta1,delta2)
	
	tl2 <- t(delta)%*%(solve(Ltt)-Ltt2)%*%delta
	tl2 <- tl2/sqrt(sum(diag(t(tl2)%*%tl2)))
	pr <- abs(eigen(tl2,symmetric=T)$vector[,length(y)])
	pr2 <- abs(diag(tl2))
	
	if(q>0){
      if(l>0){
      salida <- list(p=p,l=l,q=q,qm=0,gle.phi=gle,coefs.mu=ajuste$theta[1:p],se.mu=sqrt(diag(vcov.beta)),filas=names(beta0),coefs.phi=ajuste$theta[(p+1):(p+l+q)],se.phi=se.psi,
	  filas2=colnames(W),z.hat=z.hat,xi=xi,lambda.phi=lambda,family=family,v=v(z.hat),dg=dg,fg=fg, cdfz=cdfz(z.hat),mu.fitted=mues, phi.fitted=phi.es,deviance.mu=deviance.mu(z.hat),
	  deviance.phi=deviance.phi(z.hat),cw=cbind(cw,cw2),pr=cbind(pr,pr2),h=h,lpdf=lpdf(z.hat,phi.es),xix=xix,call="")
	  }
      if(l==0){
      salida <- list(p=p,l=l,q=q,qm=0,gle.phi=gle,coefs.mu=ajuste$theta[1:p],se.mu=sqrt(diag(vcov.beta)),filas=names(beta0),z.hat=(y-mues)/sqrt(phi.es),xi=xi,lambda.phi=lambda,
	  family=family,v=v(z.hat),dg=dg,fg=fg,cdfz=cdfz(z.hat),deviance.mu=deviance.mu(z.hat),mu.fitted=mues, phi.fitted=phi.es,deviance.phi=deviance.phi(z.hat),coefs.phi=ajuste$theta[(p+1):(p+q)],se.phi=se.psi,
	  cw=cbind(cw,cw2),pr=cbind(pr,pr2),h=h,lpdf=lpdf(z.hat,phi.es),xix=xix,call="")
	  }
    }
    else{
      if(l>0){
        salida <- list(p=p,l=l,q=q,qm=0,coefs.mu=ajuste$theta[1:p],se.mu=sqrt(diag(vcov.beta)),filas=names(beta0),coefs.phi=ajuste$theta[(p+1):(p+l)],se.phi=se.psi,
		filas2=colnames(W),z.hat=z.hat,xi=xi,family=family,v=v(z.hat),dg=dg,fg=fg,cdfz=cdfz(z.hat),mu.fitted=mues, phi.fitted=phi.es,deviance.mu=deviance.mu(z.hat),
		deviance.phi=deviance.phi(z.hat),cw=cbind(cw,cw2),pr=cbind(pr,pr2),h=h,lpdf=lpdf(z.hat,phi.es),xix=xix,call="")
	  }
    }
}
else{
	if(q>0){
	 if(l>0){
	 salida <- list(p=p,l=l,q=q,qm=0,gle.phi=gle,coefs.mu=ajuste$theta[1:p],se.mu=sqrt(diag(vcov.beta)),filas=names(beta0),coefs.phi=ajuste$theta[(p+1):(p+l+q)],se.phi=se.psi,
	 filas2=colnames(W),z.hat=z.hat,xi=xi,lambda.phi=lambda,family=family,v=v(z.hat),dg=dg,fg=fg,cdfz=cdfz(z.hat),mu.fitted=mues, phi.fitted=phi.es,deviance.mu=deviance.mu(z.hat),
	 deviance.phi=deviance.phi(z.hat),h=h,lpdf=lpdf(z.hat,phi.es),xix=xix,call="")
		}
	 if(l==0){
	 salida <- list(p=p,l=l,q=q,qm=0,gle.phi=gle,coefs.mu=ajuste$theta[1:p],se.mu=sqrt(diag(vcov.beta)),filas=names(beta0),z.hat=(y-mues)/sqrt(phi.es),xi=xi,lambda.phi=lambda,
	 family=family,v=v(z.hat),dg=dg,fg=fg,cdfz=cdfz(z.hat),deviance.mu=deviance.mu(z.hat),mu.fitted=mues, phi.fitted=phi.es,deviance.phi=deviance.phi(z.hat), coefs.phi=ajuste$theta[(p+1):(p+q)],se.phi=se.psi,
	 h=h,lpdf=lpdf(z.hat,phi.es),xix=xix,call="")
		}
	}
	else{
	 if(l>0){
	 salida <- list(p=p,l=l,q=q,qm=0,coefs.mu=ajuste$theta[1:p],se.mu=sqrt(diag(vcov.beta)),filas=names(beta0),coefs.phi=ajuste$theta[(p+1):(p+l)],se.phi=se.psi,filas2=colnames(W),
	 z.hat=z.hat,xi=xi,family=family,v=v(z.hat),dg=dg,fg=fg,cdfz=cdfz(z.hat),mu.fitted=mues, phi.fitted=phi.es,deviance.mu=deviance.mu(z.hat),deviance.phi=deviance.phi(z.hat),
     h=h,lpdf=lpdf(z.hat,phi.es),xix=xix,call="")
		}
	}
}
	
 class(salida) <- "ssym"
 salida$call <- match.call()

salida		
}
