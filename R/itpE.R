itpE <-
function(vP,objeto){
	epsilon <- objeto$epsilon
	maxiter <- objeto$maxiter	
	theta_new <- vP
	tol <- 1
	cont <- 0
	while(tol > epsilon){
		 theta <- theta_new
		 objeto$theta_work <- theta
		 theta_new <- itpM(objeto)
		 theta_pd <- ifelse(theta==rep(0,length(theta)),1,theta)
		 tol <- max(abs(theta_new-theta)/abs(theta_pd))
		 cont <- cont + 1
		 if(cont > maxiter) {
		   stop("no convergence was obtained!!",call.=FALSE) 
		 }  
    }
	list(theta=theta_new,criterion=theta,iter=cont)
}
