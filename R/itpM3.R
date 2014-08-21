itpM3 <-
function(objeto){
    p <- objeto$p
	q <- objeto$q
	qm <- objeto$qm
	l <- objeto$l
	y <- objeto$y
	v <- objeto$v
	W_bar <- objeto$W_bar
	K_psi2 <- objeto$K_psi2
	theta_new <- objeto$theta_work
	tol_EM <- 1
	cont <- 0

if(p>0){
   if(qm==0){
	 if(l>0){
	   if(q>0){
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
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
		    D2 <- D*matrix(v(z)/phi,length(phi),p)
			beta_new <- beta + (1/objeto$dg)*(1/objeto$dg)*solve(t(D)%*%D2)%*%t(D2)%*%(y-mus)
		    psi_new <- psi +  (2/(objeto$fg-1))*(2/(objeto$fg-1))*K_psi2%*%U_psi
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
  	    while(tol_EM > objeto$epsilon){
			 theta <- theta_new
	  	     beta <- theta[1:p]
    	     mus <- objeto$mu(beta)[objeto$subset]
    	     gammav <- theta[(p+1):(p+l)]
		     phi <- exp(objeto$W%*%gammav)
			 z <- (y-mus)/sqrt(phi)			 
			 psi <- theta[(p+1):length(theta)]
			 U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
    		 D <- objeto$GradD(theta)[objeto$subset,]
			 D2 <- D*matrix(v(z)/phi,length(phi),p)
			 beta_new <- beta + (1/objeto$dg)*solve(t(D)%*%D2)%*%t(D2)%*%(y-mus)
		     psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
	    while(tol_EM > objeto$epsilon){
			theta <- theta_new
			beta <- theta[1:p]
		    mus <- objeto$mu(beta)[objeto$subset]
			h <- theta[(p+1):length(theta)]
		    phi <- exp(objeto$N%*%h)
			z <- (y-mus)/sqrt(phi)
			psi <- theta[(p+1):length(theta)]
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			D <- objeto$GradD(theta)[objeto$subset,]
			D2 <- D*matrix(v(z)/phi,length(phi),p)
			beta_new <- beta + (1/objeto$dg)*solve(t(D)%*%D2)%*%t(D2)%*%(y-mus)
			psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
	}
	if(qm>0){
	D_bar <- cbind(objeto$GradD(objeto$theta),objeto$Nm)
	Kpqm <- matrix(0,p+qm,p+qm)
	Kpqm[(p+1):(p+qm),(p+1):(p+qm)] <- objeto$lambda.mu*objeto$Mm
	if(l>0){
	  if(q>0){
	    while(tol_EM > objeto$epsilon){
			theta <- theta_new
		  	beta <- theta[1:p]
		    mus <- objeto$mu(beta)
	    	gammav <- theta[(p+qm+1):(p+qm+l)]
		    h <- theta[(p+qm+l+1):length(objeto$theta)]
		    hm <- theta[(p+1):(p+qm)]			
			phi <- exp(objeto$W%*%gammav + objeto$N%*%h)
			z <- (y-mus-objeto$Nm%*%hm)/sqrt(phi)			
			psi <- theta[(p+qm+1):length(theta)]
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
		    D2 <- D_bar*matrix(v(z)/phi,length(phi),p+qm)
			beta_new <- theta[1:(p+qm)] + (1/objeto$dg)*solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D2)%*%(y-mus-objeto$Nm%*%hm) - rbind(matrix(0,p,1),objeto$lambda.mu*objeto$Mm%*%hm))
		    psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
  	    while(tol_EM > objeto$epsilon){
			 theta <- theta_new
	  	     beta <- theta[1:p]
    		 hm <- theta[(p+1):(p+qm)]			 
    	     mus <- objeto$mu(beta)
    	     gammav <- theta[(p+1+qm):(p+l+qm)]
		     phi <- exp(objeto$W%*%gammav)
			 z <- (y-mus-objeto$Nm%*%hm)/sqrt(phi)			 
			 psi <- theta[(p+qm+1):length(theta)]
			 U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			 D2 <- D_bar*matrix(v(z)/phi,length(phi),p+qm)
			 beta_new <- theta[1:(p+qm)] + (1/objeto$dg)*solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D2)%*%(y-mus-objeto$Nm%*%hm) - rbind(matrix(0,p,1),objeto$lambda.mu*objeto$Mm%*%hm))
		     psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
	    while(tol_EM > objeto$epsilon){
			theta <- theta_new
			beta <- theta[1:p]
		    mus <- objeto$mu(beta)
			h <- theta[(p+qm+1):length(theta)]
		    hm <- theta[(p+1):(p+qm)]		
		    phi <- exp(objeto$N%*%h)
			z <- (y-mus-objeto$Nm%*%hm)/sqrt(phi)
			psi <- theta[(p+qm+1):length(theta)]
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			D2 <- D_bar*matrix(v(z)/phi,length(phi),p+qm)
			beta_new <- theta[1:(p+qm)] + (1/objeto$dg)*solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D2)%*%(y-mus-objeto$Nm%*%hm) - rbind(matrix(0,p,1),objeto$lambda.mu*objeto$Mm%*%hm))
			psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
	}
}
if(p==0){
	D_bar <- objeto$Nm
	Kpqm <- matrix(0,p+qm,p+qm)
	Kpqm[(p+1):(p+qm),(p+1):(p+qm)] <- objeto$lambda.mu*objeto$Mm
	if(l>0){
	  if(q>0){
	    while(tol_EM > objeto$epsilon){
			theta <- theta_new
	    	gammav <- theta[(p+qm+1):(p+qm+l)]
		    h <- theta[(p+qm+l+1):length(objeto$theta)]
		    hm <- theta[(p+1):(p+qm)]			
			phi <- exp(objeto$W%*%gammav + objeto$N%*%h)
			z <- (y-objeto$Nm%*%hm)/sqrt(phi)			
			psi <- theta[(p+qm+1):length(theta)]
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
		    D2 <- D_bar*matrix(v(z)/phi,length(phi),p+qm)
			beta_new <- theta[1:(p+qm)] + (1/objeto$dg)*solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D2)%*%(y-objeto$Nm%*%hm) - objeto$lambda.mu*objeto$Mm%*%hm)
		    psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
  	    while(tol_EM > objeto$epsilon){
			 theta <- theta_new
    		 hm <- theta[(p+1):(p+qm)]			 
    	     gammav <- theta[(p+1+qm):(p+l+qm)]
		     phi <- exp(objeto$W%*%gammav)
			 z <- (y-objeto$Nm%*%hm)/sqrt(phi)			 
			 psi <- theta[(p+qm+1):length(theta)]
			 U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1)
			 D2 <- D_bar*matrix(v(z)/phi,length(phi),p+qm)
			 beta_new <- theta[1:(p+qm)] + (1/objeto$dg)*solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D2)%*%(y-objeto$Nm%*%hm) - objeto$lambda.mu*objeto$Mm%*%hm)
		     psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
	    while(tol_EM > objeto$epsilon){
			theta <- theta_new
			h <- theta[(p+qm+1):length(theta)]
		    hm <- theta[(p+1):(p+qm)]		
		    phi <- exp(objeto$N%*%h)
			z <- (y-objeto$Nm%*%hm)/sqrt(phi)
			psi <- theta[(p+1+qm):length(theta)]
			U_psi <- (1/2)*t(W_bar)%*%(v(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			D2 <- D_bar*matrix(v(z)/phi,length(phi),p+qm)
			beta_new <- theta[1:(p+qm)] + (1/objeto$dg)*solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D2)%*%(y-objeto$Nm%*%hm) - objeto$lambda.mu*objeto$Mm%*%hm)
			psi_new <- psi +  (2/(objeto$fg-1))*K_psi2%*%U_psi
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
}
	theta_new
}
