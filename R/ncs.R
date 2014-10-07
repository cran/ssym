ncs <-
function(xx, lambda){
	if(!is.numeric(xx)) stop("Variable in natural cubic spline must be numeric!!",call.=FALSE)

	xx <- as.matrix(xx)
	x <- as.matrix(as.numeric(levels(factor(xx))))
	if(!missingArg(lambda)){
	  	if(lambda<=0) stop("Smoothing parameter must be a positive value!!",call.=FALSE)
		else status <- "known"
	}else{status <- "unknown"
	      lambda <- 1
	}

	m <- length(x)
	n <- length(xx)
	
	h <- matrix(0,m-1,1)
	Q <- matrix(0,m,m-2)
	R <- matrix(0,m-2,m-2)
	
	for(i in 1:(m-1)){
	   h[i] <- x[i+1]-x[i]
	}
	for(j in 2:(m-1)){
	   Q[j-1,j-1] <- 1/h[j-1]
	}
	for(j in 2:(m-1)){
	   Q[j,j-1] <- -1/h[j-1] -1/h[j]
	}
	for(j in 2:(m-1)){
	   Q[j+1,j-1] <- 1/h[j]
	}
	for(j in 2:(m-1)){
	   R[j-1,j-1] <- (h[j-1]+h[j])/3
	}
	for(j in 2:(m-2)){
	   R[j-1,j] <- (h[j])/6
	}
	for(j in 2:(m-2)){
	   R[j,j-1] <- (h[j])/6
	}
	K <- Q%*%solve(R)%*%t(Q)
	N <- matrix(0,n,m)
	for(i in 1:n){
	   N[i,] <- xx[i] == x
	}
attr(xx,"Q") <- Q
attr(xx,"R") <- R
attr(xx,"K") <- K
attr(xx,"N") <- N
attr(xx,"status") <- status
attr(xx,"lambda") <- lambda
xx
}
