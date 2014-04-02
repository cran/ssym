itpM <-
function(objeto){
    p <- objeto$p
	q <- objeto$q
	l <- objeto$l
	y <- objeto$y
	v <- objeto$v
	W_bar <- objeto$W_bar
	K_psi2 <- objeto$K_psi2
	theta_new <- objeto$theta_work
	tol_EM <- 1
	cont <- 0
	if(l>0){
	  if(q>0){
	  	mus_work <- objeto$mu(objeto$theta_work[1:p])[objeto$subset]			
	    gammav_work <- objeto$theta_work[(p+1):(p+l)]			
		h_work <- objeto$theta_work[(p+l+1):length(objeto$theta_work)]
		phi_work <- exp(objeto$W%*%gammav_work + objeto$N%*%h_work)			
	    z_work <- (y-mus_work)/sqrt(phi_work)
	    while(tol_EM > objeto$epsilon){
			theta <- theta_new
		  	beta <- theta[1:p]
		    mus <- objeto$mu(beta)[objeto$subset]
	    	gammav <- theta[(p+1):(p+l)]
			D <- objeto$GradD(theta)[objeto$subset,]
		    h <- theta[(p+l+1):length(theta)]
			phi <- exp(objeto$W%*%gammav + objeto$N%*%h)
			z <- (y-mus)/sqrt(phi)			
			psi <- theta[(p+1):length(theta)]
			U_psi <- (1/2)*t(W_bar)%*%(v(z_work)*z^2-1)
			U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda*objeto$M%*%h
		    D2 <- D*matrix(v(z_work)/phi,length(phi),p)
			beta_new <- beta + solve(t(D)%*%D2)%*%t(D2)%*%(y-mus)
		    psi_new <- psi +  K_psi2%*%U_psi
    		theta_new <- c(t(beta_new),t(psi_new))
    		theta_pd <- ifelse(theta==rep(0,length(theta)),1,theta)
		    tol_EM <- max(abs(theta_new-theta)/abs(theta_pd))
			cont <- cont + 1
			if(cont > objeto$maxiter) {
			   stop("no convergence was obtained!!",call.=FALSE) 
			}  
		}
	  }
	  if(q==0){
    	mus_work <- objeto$mu(objeto$theta_work[1:p])[objeto$subset]			
	    gammav_work <- objeto$theta_work[(p+1):(p+l)]			
		phi_work <- exp(objeto$W%*%gammav_work)
		z_work <- (y-mus_work)/sqrt(phi_work)
  	    while(tol_EM > objeto$epsilon){
			 theta <- theta_new
	  	     beta <- theta[1:p]
    	     mus <- objeto$mu(beta)[objeto$subset]
    	     gammav <- theta[(p+1):(p+l)]
		     phi <- exp(objeto$W%*%gammav)
			 z <- (y-mus)/sqrt(phi)			 
			 psi <- theta[(p+1):length(theta)]
			 U_psi <- (1/2)*t(W_bar)%*%(v(z_work)*z^2-1)
    		 D <- objeto$GradD(theta)[objeto$subset,]
			 D2 <- D*matrix(v(z_work)/phi,length(phi),p)
			 beta_new <- beta + solve(t(D)%*%D2)%*%t(D2)%*%(y-mus)
		     psi_new <- psi +  K_psi2%*%U_psi
     		 theta_new <- c(t(beta_new),t(psi_new))
    		 theta_pd <- ifelse(theta==rep(0,length(theta)),1,theta)
	    	 tol_EM <- max(abs(theta_new-theta)/abs(theta_pd))
			 cont <- cont + 1
			 if(cont > objeto$maxiter) {
			   stop("no convergence was obtained!!",call.=FALSE) 
			 }  
		}
	  }
	}
	if(l==0){
	  if(q>0){
	  	mus_work <- objeto$mu(objeto$theta_work[1:p])[objeto$subset]
		h_work <- objeto$theta_work[(p+1):length(objeto$theta_work)]			
		phi_work <- exp(objeto$N%*%h_work)			
		z_work <- (y-mus_work)/sqrt(phi_work)
	    while(tol_EM > objeto$epsilon){
			theta <- theta_new
			beta <- theta[1:p]
		    mus <- objeto$mu(beta)[objeto$subset]
			h <- theta[(p+1):length(theta)]
		    phi <- exp(objeto$N%*%h)
			z <- (y-mus)/sqrt(phi)
			psi <- theta[(p+1):length(theta)]
			U_psi <- (1/2)*t(W_bar)%*%(v(z_work)*z^2-1) - objeto$lambda*objeto$M%*%h
			D <- objeto$GradD(theta)[objeto$subset,]
			D2 <- D*matrix(v(z_work)/phi,length(phi),p)
			beta_new <- beta + solve(t(D)%*%D2)%*%t(D2)%*%(y-mus)
			psi_new <- psi +  K_psi2%*%U_psi
    		theta_new <- c(t(beta_new),t(psi_new))
   		    theta_pd <- ifelse(theta==rep(0,length(theta)),1,theta)
		    tol_EM <- max(abs(theta_new-theta)/abs(theta_pd))
			cont <- cont + 1
			if(cont > objeto$maxiter) {
			  stop("no convergence was obtained!!",call.=FALSE) 
			}
		}	
	  }
    }
	theta_new
}
