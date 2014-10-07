np.graph <-
function(object, which, exp, xlab, ylab, main){

	if(missingArg(exp))exp <- FALSE
	if(missingArg(xlab) || !is.character(xlab))xlab <- " "
	if(missingArg(ylab) || !is.character(ylab))ylab <- " "
	if(missingArg(main) || !is.character(main))main <- " "


	if(missingArg(which))
	stop("The value of the which argument is missing!!",call.=FALSE)

	if(which!=1 & which!=2)
	stop("value of the which argument provided by the user is invalid!!",call.=FALSE)
	
	if(which==1){
		if(object$qm == 0) stop("There are not nonparametric component in the Median/Location submodel!!",call.=FALSE)
	xx.p <- object$np.mu
	g <- object$coefs.mu[(object$p+1):length(object$coefs.mu)]
	resp <- object$z.hat*sqrt(object$phi.fitted) + attr(xx.p,"N")%*%g
	se <- object$vcov.mu[(object$p+1):length(object$coefs.mu),(object$p+1):length(object$coefs.mu)]
	typesp <- object$type.mu
	}
	if(which==2){
		if(object$q == 0) stop("There are not nonparametric component in the Skewness/Dispersion submodel!!",call.=FALSE)
	xx.p <- object$np.phi
	g <- object$coefs.phi[(object$l+1):length(object$coefs.phi)]
	resp <- log((object$z.hat*sqrt(object$phi.fitted))^2/object$xix) - log(object$phi.fitted) + attr(xx.p,"N")%*%g
	se <- object$vcov.phi[(object$l+1):length(object$coefs.phi),(object$l+1):length(object$coefs.phi)]
	typesp <- object$type.phi
	}

	if(typesp==1){
		x.p <- as.numeric(levels(factor(xx.p)))
		se <- sqrt(diag(se))
	    cont <- 1000
		ss <- ncs(xx.p)
	    gam <- solve(attr(xx.p,"R"))%*%t(attr(xx.p,"Q"))%*%g
		gam <- rbind(0,gam,0)
		t <- seq(min(x.p),max(x.p),length=cont)
		nodo <- 1
		spl <- matrix(0,cont,1)
	
		for(i in 1:cont){
			if(t[i] > x.p[nodo + 1]){
			  nodo <- nodo + 1
			}
			h <- x.p[nodo+1] - x.p[nodo]
			spl[i] <- ((t[i]-x.p[nodo])*g[nodo+1] + (x.p[nodo+1]-t[i])*g[nodo])/h -
			          (t[i]-x.p[nodo])*(x.p[nodo+1]-t[i])*((1 + (t[i]-x.p[nodo])/h)*gam[nodo+1]
					  + (1 + (x.p[nodo+1]-t[i])/h)*gam[nodo])/6
		}
	 }else{
		  iknots <- attr(attr(xx.p,"N"),"knots")
		  x.p <- seq(min(xx.p),max(xx.p),length=500)
		  Bs <- bs(x.p,knots=iknots,degree=3,intercept=TRUE)
		  spl <- Bs%*%g
		  g <- spl
		  se <- sqrt(diag(Bs%*%se%*%t(Bs)))
		  t <- x.p
		  }

	  if(exp==TRUE){
	  	  lims <- cbind(exp(g - 2*se),exp(g + 2*se))
	  	  plot(t, exp(spl), xlim=range(xx.p), ylim=range(exp(resp)), type="l", col="blue", xlab="", ylab="", main="")
		  par(new=TRUE)
		  plot(x.p, lims[,1], xlim=range(xx.p), ylim=range(exp(resp)), type="l", col="blue", xlab="", ylab="", main="", lty=3)
		  par(new=TRUE)
		  plot(x.p, lims[,2], xlim=range(xx.p), ylim=range(exp(resp)), type="l", col="blue", xlab="", ylab="", main="", lty=3)
		  par(new=TRUE)
          plot(xx.p, exp(resp), xlim=range(xx.p), ylim=range(exp(resp)), type="p", cex=0.3, lwd=3, xlab=xlab, ylab=ylab, main=main)

	  }else{
	  	  lims <- cbind(g - 2*se,g + 2*se)
	  	  plot(t, spl, xlim=range(xx.p), ylim=range(resp), type="l", col="blue", xlab="", ylab="", main="")
		  par(new=TRUE)
		  plot(x.p, lims[,1], xlim=range(xx.p), ylim=range(resp), type="l", col="blue", xlab="", ylab="", main="", lty=3)
		  par(new=TRUE)
		  plot(x.p, lims[,2], xlim=range(xx.p), ylim=range(resp), type="l", col="blue", xlab="", ylab="", main="", lty=3)
		  par(new=TRUE)
          plot(xx.p, resp, xlim=range(xx.p), ylim=range(resp), type="p", cex=0.3, lwd=3, xlab=xlab, ylab=ylab, main=main)
	  }
}
