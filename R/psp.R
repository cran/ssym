psp <-
function(xx, lambda, nknots, diff){
   if(!is.numeric(xx)) stop("Variable in P-spline must be numeric!!",call.=FALSE)
	if(!missingArg(lambda)){
	  	if(lambda<=0) stop("Smoothing parameter must be a positive value!!",call.=FALSE)
		else status <- "known"
	}else{status <- "unknown"
	      lambda <- 1
	}
	
   x <- as.matrix(xx)

   if(missingArg(nknots)){
	   n <- length(x)
	   dknots <- floor(n^(1/3))
   }else{if(floor(nknots)<=0) stop("Number of internal knots must be a positive integer!!",call.=FALSE)
   		 dknots <- floor(nknots)
   }
   iknots <- quantile(xx,prob=seq(0,1,length=(dknots+2)))[2:(dknots+1)]

   if(!missingArg(diff)){ 
	  	if(floor(diff)<=0) stop("Order of the difference penalty term must be a positive integer!!",call.=FALSE)
		diff <- floor(diff)
   }	
   else diff <- 2

   xx2 <- bs(xx,knots=iknots,degree=3,intercept=TRUE)
   N <- as.matrix(xx2)
   P <- diff(diag(ncol(N)),differences=diff)
   K <- t(P)%*%P

   if(sum(is.na(N))>0) stop("Too many internal knots!!",call.=FALSE)
   
attr(x,"K") <- K
attr(x,"N") <- N
attr(x,"status") <- status
attr(x,"lambda") <- lambda
x
}
