ssym.nl <-
function(formula, start, family, xi, data, epsilon, maxiter, subset, link.phi, local.influence, spec, std.out){

if(family!="Normal" & family!="Slash" & family!="Hyperbolic" & family!="Sinh-t" & family!="Sinh-normal" & family!="Contnormal" & family!="Powerexp" & family!="Student")
stop("family of distributions specified by the user is not supported!!",call.=FALSE)

if(missingArg(link.phi)) link.phi <- "log"
if(link.phi=="logarithmic") link.phi <- "log"
if(link.phi!="identity" & link.phi!="log") 
stop("link.phi function specified by the user is not supported!!",call.=FALSE)
if(missingArg(std.out) || std.out!="FALSE") std.out <- "TRUE"


if(family=="Slash" | family=="Hyperbolic" | family=="Sinh-t" | family=="Sinh-normal" | family=="Contnormal" | family=="Powerexp" | family=="Student"){
  if(missingArg(xi)) stop("for the family of distributions specified by the user an extra parameter is required!!", call.=FALSE) 
}
	 reem <- function(aa,b){
		 ag <- aa
		 ag <- sub("(", "", ag,fixed=TRUE)
		 ag <- sub(")", "", ag,fixed=TRUE)
		 ag <- sub(b, "", ag,fixed=TRUE)
		 ag <- strsplit(ag, ",")
		 ag <- ag[[1]][1]
		 ag
	 }

	 bdiag <- function(x1,x2){
		 x1 <- as.matrix(x1)
		 x2 <- as.matrix(x2)
		 dx1 <-	dim(x1)
		 dx2 <-	dim(x2)
		 x <- matrix(0,dx1[1]+dx2[1],dx1[2]+dx2[2])
		 x[1:dx1[1],1:dx1[2]] <- x1
		 x[(dx1[1]+1):(dx1[1]+dx2[1]),(dx1[2]+1):(dx1[2]+dx2[2])] <- x2		 
    x
    }
	Klambda <- function(eta,M,q){
			  pos <- c(0,cumsum(q))
			  KM <- M
		      for(i in 1:length(q)){
			     KM[(pos[i]+1):pos[i+1],(pos[i]+1):pos[i+1]] <- KM[(pos[i]+1):pos[i+1],(pos[i]+1):pos[i+1]]*exp(eta[i])
		      }
			  KM
	}

    if(missingArg(epsilon)) epsilon <- 0.0000001
    if(missingArg(maxiter)) maxiter <- 1000

	 mf <- match.call(expand.dots = FALSE)
	 m <- match(c("formula"), names(mf), 0)
	 mf <- mf[c(1, m)]

     if (missingArg(data)) 
        data <- environment(formula)

	if(sum(is.na(match(names(start),all.vars(formula))))>0)
	stop("Some of the specified parameters in the argument start are not involved in the nonlinear function to be fitted!!",call.=FALSE)

	formula <- Formula(formula)
	if (length(formula)[2L] < 2L)
        formula <- as.Formula(formula(formula), ~1)

	nl <- formula(formula, lhs=0, rhs=1)
	ll <- formula(formula, lhs=1, rhs=2)

	response <- as.matrix(eval(ll[[2]], data))
	if(missingArg(subset)) subset <- 1:length(response)
	y <- as.matrix(response[subset])
	n <- length(y)
	if(missingArg(spec)) spe <- 2
	else{if(spec!="AIC" & spec!="BIC") stop("The specified criterion to estimate the smoothing parameter is not supported!!",call.=FALSE)
		 spe <- ifelse(spec=="AIC",2,log(n))
	}
	pars <- names(start)
	p <- length(start)
	qm <- 0
    mu <- function(beta){
	    for(i in 1:p) assign(pars[i], beta[i], environment(formula))
		temp <- as.matrix(eval(nl[[2]], data))
		if(nrow(temp)>1) as.matrix(temp[subset,])
		else as.matrix(temp)
	}

    GradD <- function(vP){jacobian(mu,vP[1:p])}
    if(length(pars)==1 & nrow(mu(start))==1){
    	mu <- function(beta) matrix(beta[1],length(y),1)
    	GradD <- function(vP) matrix(1,length(y),1)
    }   
	penm <- matrix(0,p,p)
	pspm <- GradD(start)

	
	 w <- model.matrix(ll, data)
	 wa  <- "ncs("	 
	 wa2 <- "psp("	 
	 wb <- colnames(w)
	 idw  <- grepl(wa,wb,fixed=TRUE)
	 idw2 <- grepl(wa2,wb,fixed=TRUE)
	 idw3 <- idw2 + idw
	 which.phi <- ifelse(idw3 > 0,"ncs","")
	 which.phi[idw2==1] <- "psp"

	 type.phi <- matrix(1,sum(idw3),1)
	 type.phi[idw2==1] <- 2
	 pspp <- matrix(0,n,1)
	 model.matrix.phi <- matrix(0,n,1)	 
	 q <- 0
	 penp <- matrix(0,1,1)
	 lambdas.phi <- matrix(0,1,1)		
	 statusp <- "known"

	 if(sum(idw3) < length(idw3)){
		W <- as.matrix(w[subset,!idw3])
		colnames(W) <- colnames(w)[!idw3]
		l <- ncol(W)
		haver <- apply(W,2,var)
		haver2 <- colnames(W)
		if(mean(W[,1])==1) haver[1] <- 1 
		W <- as.matrix(W[,haver!=0])
		colnames(W) <- haver2[haver!=0]
		l <- ncol(W)
		pspp <- cbind(pspp,W)
		penp <- bdiag(penp,matrix(0,l,l))
		model.matrix.phi <- cbind(model.matrix.phi,W)
	 }else{l <- 0}
	 if(sum(idw3) >= 1){
		for(i in 1:length(idw3)){
			if(idw3[i]==1){
				temp <- eval(parse(text=reem(wb[i],which.phi[i])), data)
				temp <- as.matrix(temp[subset])
				xx.p <- eval(parse(text=sub(reem(wb[i],which.phi[i]),"temp",wb[i],fixed=TRUE)))
				cent <- qr.Q(qr(matrix(1,ncol(attr(xx.p,"N")),1)),complete=TRUE)[,-1]
				pspp <- cbind(pspp,attr(xx.p,"N")%*%cent)
				q <- cbind(q,ncol(attr(xx.p,"N"))-1)
				penp <- bdiag(penp,t(cent)%*%attr(xx.p,"K")%*%cent)
				if(attr(xx.p,"status")=="known") lambdas.phi <- cbind(lambdas.phi,attr(xx.p,"lambda"))
				else lambdas.phi <- cbind(lambdas.phi,0)
				statusp <- cbind(statusp,attr(xx.p,"status"))
			}	
		}
	lambdas.phi <- lambdas.phi[-1]
	statusp <- statusp[-1]
	vnpsp <- as.matrix(w[subset,idw3==1])
	colnames(vnpsp) <- colnames(w)[idw3==1]
    model.matrix.phi <- cbind(model.matrix.phi,vnpsp)	
	q <- q[-1]
	}	
	pspp <- as.matrix(pspp[,-1])
	penp <- as.matrix(penp[-1,-1])
    pspp2 <- crossprod(pspp)
	col.m2 <- colnames(model.matrix.phi)[-1]
	model.matrix.phi <- as.matrix(model.matrix.phi[,-1])
	colnames(model.matrix.phi) <- col.m2

	
if(family=="Normal"){
	xi <- 0
	v <- function(z) rep(1,length(z)) 
	vp <- function(z) rep(0,length(z))
	dg <- 1
    fg <- 3
	deviance.mu <- function(z) z^2
	deviance.phi <- function(z) abs(z^2-1-log(z^2))
	cdf <- function(z) pnorm(z)
	pdf <- function(z) dnorm(z)
	xix <- 1
}
if(family=="Student"){
	if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	nu <- xi[1]
	v <- function(z) (nu + 1)/(nu + z^2) 
	vp <- function(z) -2*z*(nu + 1)/((nu + z^2)^2)
	dg <- (nu + 1)/(nu + 3)
    fg <- 3*(nu + 1)/(nu + 3)
	deviance.mu <- function(z) abs((nu+1)*log(1 + z^2/nu))
	deviance.phi <- function(z) abs((nu+1)*log((nu + z^2)/(nu + 1)) -log(z^2))
	cdf <- function(z) pt(z,nu)
	pdf <- function(z) dt(z,nu)
	xix <- nu/(nu-2)
	if(nu<=2) xix <- as.null(xix)
}
if(family=="Contnormal"){
	if(xi[1]<=0  | xi[1]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
	if(xi[2]<=0  | xi[2]>=1) stop("the extra parameters must be within the interval (0, 1)!!",call.=FALSE)
	eta <- xi[1:2]
	v <- function(z) (eta[2]^(3/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))/(eta[2]^(1/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))
	vp <- function(z) eta[2]^(1/2)*eta[1]*exp(z^2*(1-eta[2])/2)*(1-eta[2])*z*(eta[2]*(1-eta[1]) - (1-eta[1]))/((eta[2]^(1/2)*eta[1]*exp(z^2*(1-eta[2])/2) + (1-eta[1]))^2)
	deviance.mu <- function(z) abs(-2*log(((eta[2]^(1/2)*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))/(eta[2]^(1/2)*eta[1]*dnorm(0) + (1-eta[1])*dnorm(0)))))
	tau <-  uniroot(function(x) v(x)*x^2 -1, lower=0, upper=35)$root
	deviance.phi <- function(z) abs(-2*log((abs(z/tau))*((eta[2]^(1/2)*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))/(eta[2]^(1/2)*eta[1]*dnorm(tau*sqrt(eta[2])) + (1-eta[1])*dnorm(tau)))))
	cdf <- function(z) eta[1]*pnorm(z*sqrt(eta[2])) + (1-eta[1])*pnorm(z)
	pdf <- function(z) sqrt(eta[2])*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z)
	xix <- eta[1]/eta[2] + (1-eta[1])
 	dgd <- function(z) (v(z)*z)^2*(sqrt(eta[2])*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))
	dg <- 2*integrate(dgd,0,35)$value
 	fgd <- function(z) (v(z))^2*z^4*(sqrt(eta[2])*eta[1]*dnorm(z*sqrt(eta[2])) + (1-eta[1])*dnorm(z))
	fg <- 2*integrate(fgd,0,35)$value
}

if(family=="Powerexp"){
	if(xi[1]<=-1  | xi[1]>=1) stop("the extra parameter must be within the interval (-1, 1)!!",call.=FALSE)
	kk <- xi[1]
	v <- function(z) abs(z)^(-(2*kk/(1+kk)))/(1+kk)
	vp <- function(z) -2*kk*ifelse(z>=0,1,-1)*abs(z)^(-((3*kk + 1)/(1+kk)))/(1+kk)^2
	dg <- 2^(1-kk)*gamma((3-kk)/2)/((1+kk)^2*gamma((1+kk)/2))
	fg <- (kk + 3)/(kk + 1)
	deviance.mu <- function(z) abs((abs(z))^(2/(1+kk)))
	deviance.phi <- function(z) abs((abs(z))^(2/(1+kk)) - (1+kk) - log(z^2/((1+kk)^(1+kk))))
	pp <- 2/(kk+1)
	sigmap <- (1+kk)^((kk+1)/2)
	cdf <- function(z) pnormp(z,mu=0,sigmap=sigmap,p=pp)
	pdf <- function(z) dnormp(z,mu=0,sigmap=sigmap,p=pp)
	xix <- 2^(1+kk)*gamma(3*(1+kk)/2)/gamma((1+kk)/2)
}
if(family=="Sinh-normal"){
   	if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	alpha <- xi[1]
	v <- function(z) 4*sinh(z)*cosh(z)/(alpha^2*z) - tanh(z)/z
	vp <- function(z) ((cosh(z)*z-sinh(z))/z^2)*(4*cosh(z)/alpha^2 - 1/cosh(z)) + (sinh(z)/z)*(4*sinh(z)/alpha^2 + sinh(z)/(cosh(z)^2))
	dg <- 2 + 4/(alpha^2) - (sqrt(2*pi)/alpha)*(1-2*(pnorm(sqrt(2)/alpha,mean=0,sd=sqrt(2)/2)-0.5))*exp(2/(alpha^2))
	dshn <- function(z) 2*cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)/(alpha*sqrt(2*pi))
 	fgf <- function(z) dshn(z)*(4*sinh(z)*cosh(z)/(alpha^2) - tanh(z))^2*z^2
	fg <- 2*integrate(fgf,0,60)$value
    deviance.mu <- function(z){if(alpha<=2) abs(4*(sinh(z))^2/alpha^2 - log((cosh(z))^2))
				               else{2*log(dshn(acosh(alpha/2))/dshn(z))}}
	tau <- uniroot(function(x) 4*x*sinh(x)*cosh(x)/(alpha^2) - tanh(x)*x -1, lower=0, upper=50)$root
    deviance.phi <- function(z){
	            a <- 2*log(cosh(tau)*exp(-(2/alpha^2)*(sinh(tau))^2)) + log(tau^2)                                       
				b <- 2*log(cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)) + log(z^2)
				ifelse(a<b,0,a-b)}
	cdf <- function(z) pnorm(2*sinh(z)/alpha)
	pdf <- function(z) (2*cosh(sqrt(z^2))*exp(-2*sinh(sqrt(z^2))*sinh(sqrt(z^2))/alpha^2)/(sqrt(2*pi)*alpha))
	fgf <- function(z) dshn(z)*z^2
	xix <- 2*integrate(fgf,0,20)$value
}	

if(family=="Sinh-t"){
   	if(xi[1]<=0 | xi[2]<=0) stop("the extra parameters must be positive!!",call.=FALSE)
	alpha <- xi[1]
	nu <- xi[2]	
	v <- function(z)  4*sinh(z)*cosh(z)*(nu+1)/(alpha^2*z*nu + 4*z*sinh(z)*sinh(z)) - tanh(z)/z
	vp <- function(z){
	  (4*(nu+1)/z)*((alpha^2*nu + 4*sinh(z)*sinh(z))*(sinh(z)*sinh(z) + cosh(z)*cosh(z)) - 8*sinh(z)*sinh(z)*cosh(z)*cosh(z))/(alpha^2*nu + 4*sinh(z)*sinh(z))^2  - 1/(z*cosh(z)*cosh(z))  - v(z)*z/z^2
	}
	dsht <- function(z) (gamma((nu+1)/2)/(gamma(nu/2)*sqrt(nu*pi)))*(2/alpha)*cosh(z)*(1 + 4*(sinh(z))^2/(nu*alpha^2))^(-(nu+1)/2)
 	dgd <- function(z) dsht(z)*(4*sinh(z)*cosh(z)*(nu+1)/(alpha^2*nu + 4*sinh(z)*sinh(z)) - tanh(z))^2
	dg <- 2*integrate(dgd,0,60)$value
 	fgf <- function(z) dsht(z)*(4*sinh(z)*cosh(z)*(nu+1)/(alpha^2*nu + 4*sinh(z)*sinh(z)) - tanh(z))^2*z^2
	fg <- 2*integrate(fgf,0,70)$value
    deviance.mu <- function(z){if(alpha<=2*sqrt(1 + 1/nu)) (nu + 1)*log(1 + 4*(sinh(z))^2/(nu*alpha^2)) - log((cosh(z))^2)
				               else{2*log(dsht(acosh(sqrt(alpha^2/4 - 1/nu)))/dsht(z))}}
	tau <- uniroot(function(x) 4*sinh(x)*cosh(x)*(nu+1)*x/(alpha^2*nu + 4*sinh(x)*sinh(x)) - tanh(x)*x -1, lower=0, upper=80)$root
    deviance.phi <- function(z){
	            a <- 2*log(cosh(tau)*(1 + 4*(sinh(tau))^2/(nu*alpha^2))^(-(nu+1)/2)) + log(tau^2)                                       
				b <- 2*log(cosh(z)*(1 + 4*(sinh(z))^2/(nu*alpha^2))^(-(nu+1)/2)) + log(z^2)                                       
				ifelse(a<b,0,a-b)}
	cdf <- function(z) pt(2*sinh(z)/alpha,nu)
	pdf <- function(z) dsht(z)
	fgf <- function(z) dsht(z)*z^2
	xix <- 2*integrate(fgf,0,50)$value
	dshn <- function(z) 2*cosh(z)*exp(-(2/alpha^2)*(sinh(z))^2)/(alpha*sqrt(2*pi))
 	fgf <- function(z) dshn(z)*(4*sinh(z)*cosh(z)/(alpha^2) - tanh(z))^2*z^2
}

if(family=="Hyperbolic"){
    if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	nu <- xi[1]
	v <- function(z) nu/sqrt(1 + z^2)
	vp <- function(z) -nu*z/((1 + z^2)^(3/2))
	dh <- function(z) exp(-nu*sqrt(1+z^2))/(2*besselK(nu, 1))
 	dgd <- function(z) dh(z)*(nu*z/sqrt(1 + z^2))^2
	dg <- 2*integrate(dgd,0,Inf)$value
 	fgf <- function(z) dh(z)*(nu*z^2/sqrt(1 + z^2))^2
	fg <- 2*integrate(fgf,0,Inf)$value
    deviance.mu <- function(z) abs(2*nu*(sqrt(1+z^2)-1))
	tau <- sqrt((1 + sqrt(1 + 4*nu^2))/(2*nu^2))
    deviance.phi <- function(z) abs(2*nu*(sqrt(1+z^2) - sqrt(1 + tau^2)) - log(z^2/tau^2))
	cdf <- function(z){temporal <- matrix(0,length(z),1)
	                   for(gg in 1:length(z)) temporal[gg] <- integrate(dh,-Inf,z[gg])$value
					   temporal}
	pdf <- function(z) dh(z)
	fgf <- function(z) dh(z)*z^2
	xix <- 2*integrate(fgf,0,Inf)$value
}

if(family=="Slash"){
   	if(xi[1]<=0) stop("the extra parameter must be positive!!",call.=FALSE)
	nu <- xi[1]
	G <- function(a,x) 	gamma(a)*pgamma(1,shape=a,scale=1/x)/(x^a)
	v <- function(z) G(nu+3/2,z^2/2)/G(nu+1/2,z^2/2)
	vp <- function(z) grad(v,z)
	ds <- function(z) nu*G(nu+1/2,z^2/2)/sqrt(2*pi)
	gdg <- function(z) ds(z)*(v(z))^2*z^2
	dg <- 2*integrate(gdg,0,Inf)$value
	gfg <- function(z) ds(z)*(v(z))^2*z^4
	fg <- 2*integrate(gfg,0,Inf)$value
    deviance.mu <- function(z) abs(2*log(2/(2*nu+1))-2*log(G(nu+1/2,z^2/2)))
	tau <- uniroot(function(x) v(x)*x^2 -1, lower=0.0001, upper=1000)$root
    deviance.phi <- function(z){
	            a <- 2*log(G(nu+1/2,tau^2/2)) + log(tau^2)                                       
				b <- 2*log(G(nu+1/2,z^2/2)) + log(z^2)                                                                              
				ifelse(a<b,0,a-b)}
	cdf <- function(z){temporal <- matrix(0,length(z),1)
	                   for(gg in 1:length(z)) temporal[gg] <- integrate(ds,-Inf,z[gg])$value
					   temporal}
	pdf <- function(z) ds(z)
	gfg <- function(z){ds(z)*z^2}
	xix <- 2*integrate(gfg,0,60)$value
}

l.mu <- function(y) y
l.mu.i <- function(y) y
l1.mu <- function(mu) matrix(1,length(mu),1)
attr(l1.mu,"link") <- "identity"

if(missingArg(link.phi) || link.phi=="log"){
	l.phi <- function(y) log(y)
	l.phi.i <- function(y) exp(y)
	l1.phi <- function(phi) matrix(1,length(phi),1)
	l2.phi <- function(phi) -1/phi^2	
	attr(l1.phi,"link") <- "logarithmic"}
else{
	l.phi <- function(y) y
	l.phi.i <- function(y) y
	l1.phi <- function(phi) 1/phi
	l2.phi <- function(phi) matrix(0,length(phi),1)	
	attr(l1.phi,"link") <- "identity"}
	
objeto <- list(maxiter=maxiter, l.mu=l.mu, l.phi=l.phi, l.mu.i=l.mu.i, l.phi.i=l.phi.i, l1.mu=l1.mu, l1.phi=l1.phi, mu=mu, GradD=GradD, epsilon=epsilon, y=y, family=family, xi=xi, dg=dg, fg=fg, v=v, l=l, q=q, p=p, qm=qm, pspp=pspp, pspm=pspm, penm=penm, orig="nonlinear")
thetam0 <- start
thetap0 <- solve(crossprod(pspp) + penp)%*%crossprod(pspp,l.phi((y-mu(thetam0))^2/ifelse(is.null(xix),1,xix)))
theta0 <- c(thetam0,thetap0)

if(any(statusp == "unknown")){
	frlambda <- function(lambdas){
                lambdas.phi[statusp == "unknown"]	<- exp(lambdas[1:sum(statusp == "unknown")])
				if(l>0) penp2 <- Klambda(c(0,log(lambdas.phi)),penp,c(l,q))
				else  penp2 <- Klambda(log(lambdas.phi),penp,q)

				objeto$penp <- penp2

				if(family=="Sinh-t"){
				vP <- itpE2(theta0,objeto)}
				else{if( (family=="Powerexp" && xi < 0) || (family=="Sinh-normal") || (family=="Normal") ) vP <- itpE3(theta0,objeto)
				     else vP <- itpE(theta0,objeto)
				}
				
				mu_es <- mu(vP[1:(p+sum(qm))])
				phi_es <- l.phi.i(pspp%*%vP[(p+sum(qm)+1):(p+sum(qm)+l+sum(q))])
				z_es <- (y - mu_es)/sqrt(phi_es)

			    if(attr(l1.phi,"link")!="logarithmic") pspp2 <- crossprod(pspp*matrix(l1.phi(phi_es),n,ncol(pspp)))
			    gle <- sum(diag(solve(pspp2 + (4/(fg-1))*penp2)%*%pspp2))
				
				-2*sum(log(pdf(z_es)) -(1/2)*log(phi_es)) + spe*gle
	}
	  salida <- try(nlminb(rep(0,sum(statusp == "unknown")),frlambda), silent=TRUE)
	  if(!is.list(salida)){
	    salida <- try(nlm(frlambda,rep(0,sum(statusp == "unknown"))), silent=TRUE)
		if(!is.list(salida)){
		  salida <- optim(rep(0,sum(statusp == "unknown")),frlambda,method="BFGS")
		  par <- salida$par
		}else{par <- salida$estimate}
	  }else{par <- salida$par}

	  
	  if(any(statusp == "unknown")) lambdas.phi <- exp(par[1:sum(statusp == "unknown")])
}	

	if(sum(q)>0){
	  if(l>0) penp <- Klambda(c(0,log(lambdas.phi)),penp,c(l,q))
	  else  penp <- Klambda(log(lambdas.phi),penp,q)
	}
	objeto$penm <- penm
	objeto$penp <- penp

	
if(family=="Sinh-t"){
vP <- itpE2(theta0,objeto)}
else{if( (family=="Powerexp" && xi < 0) || (family=="Sinh-normal") || (family=="Normal") ) vP <- itpE3(theta0,objeto)
     else vP <- itpE(theta0,objeto)
}

thetam <- vP[1:p]
thetap <- vP[(p+1):(p+l+sum(q))]
mu_es <- mu(thetam)
phi_es <- l.phi.i(pspp%*%thetap)
z_es <- (y - mu_es)/sqrt(phi_es)
v_es <- v(z_es)

if(std.out=="TRUE"){
	pspm <- GradD(thetam)
	score.mu <- crossprod(pspm,z_es*v_es/sqrt(phi_es))
	score.phi <- crossprod(pspp,l1.phi(phi_es)*(v_es*z_es^2 - 1)/2) - penp%*%thetap
	vcov.mu <- solve(dg*crossprod(pspm,pspm*matrix(l1.mu(z_es)^2/phi_es,n,ncol(pspm))))
	if(attr(l1.phi,"link")!="logarithmic") pspp2 <- crossprod(pspp,pspp*matrix(l1.phi(phi_es)^2,n,l+sum(q)))
	vcov.phi <- solve(((fg-1)/4)*pspp2 + penp)
	if(sum(q)>0) vcov.phi <- crossprod(vcov.phi,tcrossprod(((fg-1)/4)*pspp2,vcov.phi))
	
	gle.mu <- p
	if(sum(q) > 0 ){
	  stes.phi <- matrix(-1,length(q),2)
	  pos <- c(0,cumsum(q)) + l
      hatm <- diag(tcrossprod(solve(pspp2 + (4/(fg-1))*penp),pspp2))
	  if(l>0){
	    gle.phi <- matrix(0,length(q)+1,1)
	  	gle.phi[1] <- sum(hatm[1:l])
	  }
	  else gle.phi <- matrix(0,length(q),1)
	  for(i in 1:length(q)){
		  gle.phi[i+min(1,l)] <- sum(hatm[(pos[i]+1):(pos[i+1])])
		  invm <- try(solve(vcov.phi[(pos[i]+1):pos[i+1],(pos[i]+1):pos[i+1]]), silent=TRUE)
		  if(!is.matrix(invm)) invm <- solve(vcov.phi[(pos[i]+1):pos[i+1],(pos[i]+1):pos[i+1]] + diag(q[i])*0.000000001)
		  stes.phi[i,1] <- crossprod(thetap[(pos[i]+1):pos[i+1]],invm%*%thetap[(pos[i]+1):pos[i+1]])
		  stes.phi[i,2] <- 1 - pchisq(stes.phi[i,1],df=q[i])
	  }
	}else gle.phi <- l
	  
	out_ <- list(p=p,qm=qm,l=l,q=q,y=y,score.mu=score.mu,score.phi=score.phi,vcov.mu=vcov.mu,vcov.phi=vcov.phi,family=family,xi=xi,mu=mu, GradD=GradD,
	               model.matrix.mu=GradD(thetam), model.matrix.phi=model.matrix.phi, deviance.mu=deviance.mu(z_es), deviance.phi=deviance.phi(z_es), deviance.mu.f=deviance.mu, deviance.phi.f=deviance.phi,
				   z_es=z_es, theta.mu=thetam, theta.phi=thetap, cdfz=cdf(z_es), lpdf=(log(pdf(z_es))-(1/2)*log(phi_es)), xix=xix, weights=v_es/dg, censored="FALSE",
				   dg=dg, fg=fg, mu.fitted=mu_es, phi.fitted=phi_es, v=v, orig="nonlinear", pspm=pspm, pspp=pspp, penm=penm, penp=penp, gle.mu=gle.mu, gle.phi=gle.phi,
				   l.mu=l.mu, l.phi=l.phi, l.mu.i=l.mu.i, l.phi.i=l.phi.i, l1.mu=l1.mu, l1.phi=l1.phi)
	
	colnames(out_$model.matrix.mu) <- names(start)			   
	if(sum(q)>0){
	  out_$which.phi <- which.phi
	  out_$lambdas.phi <- lambdas.phi
	  out_$stes.phi <- stes.phi
	  out_$type.phi <- type.phi
	}

	if(!missingArg(local.influence)){
	   splb <- matrix(0,p,p)
       for(ij in 1:length(phi_es)){
   		   mu_ij <- function(vP) mu(vP)[ij]
   		   splb <- splb + (z_es[ij]*v(z_es[ij])/sqrt(phi_es[ij]))*hessian(func=mu_ij,x=thetam)
        }
		vp_es <- vp(z_es)
		Dc   <- (vp_es*z_es + v_es)
		Dcar <- (vp_es*z_es^2 + 2*z_es*v_es)/2
		Dcab <- (Dcar*z_es + (v_es*z_es^2 - 1)*(1 + phi_es^2*l2.phi(phi_es)*l1.phi(phi_es)))/2
		
		Ltt <- matrix(0,ncol(pspm)+ncol(pspp),ncol(pspm)+ncol(pspp))
		Ltt2 <- matrix(0,ncol(pspm)+ncol(pspp),ncol(pspm)+ncol(pspp))	
		Ltt[1:ncol(pspm),1:ncol(pspm)] <- -crossprod(pspm,pspm*matrix(Dc/phi_es,n,ncol(pspm))) + splb
		Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)] <- -crossprod(pspp,pspm*matrix(l1.phi(phi_es)*Dcar/sqrt(phi_es),n,ncol(pspm)))
		Ltt[1:ncol(pspm),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- t(Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)])
		Lgg <- -crossprod(pspp,pspp*matrix(l1.phi(phi_es)^2*Dcab,n,ncol(pspp))) - penp
		Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- Lgg
		Ltt2[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- solve(Lgg)
	
		delta <- rbind(t(pspm*matrix(v_es*z_es/sqrt(phi_es),length(phi_es),ncol(pspm))), (1/2)*t(matrix(l1.phi(phi_es)*(v_es*z_es^2-1),length(phi_es),ncol(pspp))*pspp))
		tl <- t(delta)%*%(solve(Ltt)-Ltt2)%*%delta
		tl <- tl/sqrt(sum(diag(t(tl)%*%tl)))
		cw <- abs(eigen(tl,symmetric=TRUE)$vector[,length(phi_es)])
		cw2 <- abs(diag(tl))
		tl <- t(delta)%*%(solve(Ltt))%*%delta
		tl <- tl/sqrt(sum(diag(t(tl)%*%tl)))
		cw.theta <- abs(eigen(tl,symmetric=TRUE)$vector[,length(y)])
		cw2.theta <- abs(diag(tl))
		
	    delta <- rbind(t(pspm*matrix(Dc/phi_es,length(phi_es),ncol(pspm))),t(matrix(l1.phi(phi_es)*Dcar/sqrt(phi_es),length(phi_es),ncol(pspp))*pspp))
		tl2 <- t(delta)%*%(solve(Ltt)-Ltt2)%*%delta
		tl2 <- tl2/sqrt(sum(diag(t(tl2)%*%tl2)))
		pr <- abs(eigen(tl2,symmetric=TRUE)$vector[,length(y)])
		pr2 <- abs(diag(tl2))
		tl2 <- t(delta)%*%(solve(Ltt))%*%delta
		tl2 <- tl2/sqrt(sum(diag(t(tl2)%*%tl2)))
		pr.theta <- abs(eigen(tl2,symmetric=TRUE)$vector[,length(y)])
		pr2.theta <- abs(diag(tl2))
	
		out_$cw <- cbind(cw,cw2)
		out_$pr <- cbind(pr,pr2)
		out_$cw.theta <- cbind(cw.theta,cw2.theta)
		out_$pr.theta <- cbind(pr.theta,pr2.theta)
}
}else{out_ <- list(cdfz=cdf(z_es), lpdf=(log(pdf(z_es))-(1/2)*log(phi_es)), censored="FALSE")}

 class(out_) <- "ssym"
 out_$call <- match.call()

 out_  

}
