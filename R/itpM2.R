itpM2 <-
function(objeto){
    p <- objeto$p
    qm <- objeto$qm
	q <- objeto$q
	l <- objeto$l
	y <- objeto$y
	v <- objeto$v
	alpha <- objeto$xi[1]
	dg2 <- function(alphas){
			2 + 4/(alphas^2) - (sqrt(2*pi)/alphas)*(1-2*(pnorm(sqrt(2)/alphas,mean=0,sd=sqrt(2)/2)-0.5))*exp(2/(alphas^2))
	}
	fg2 <- function(alphas){
		  (2.071 + alphas*(1.936)      + alphas^2*(-6.003e-02) + alphas^3*(1.337e-03)  + alphas^4*(-1.748e-05) + alphas^5*(1.202e-07) + alphas^6*(-3.343e-10))*(alphas>2) + 
		  (3.001 + alphas*(-1.030e-02) + alphas^2*(1.050)      + alphas^3*(-8.828e-02) + alphas^4*(-1.859e-01) + alphas^5*(8.899e-02) + alphas^6*(-1.327e-02))*(alphas<=2)
	}
	family <- objeto$family
	W_bar <- objeto$W_bar
	theta_new <- objeto$theta_work
	tol_EM <- 1
	cont <- 0

if(p>0){
	if(qm==0){
		if(l>0){
		  if(q>0){
			if(family=="Sinh-t"){
				mus_work <- objeto$mu(objeto$theta_work[1:p])[objeto$subset]			
			    gammav_work <- objeto$theta_work[(p+1):(p+l)]			
				h_work <- objeto$theta_work[(p+l+1):length(objeto$theta_work)]
				phi_work <- exp(objeto$W%*%gammav_work + objeto$N%*%h_work)			
			    z_work <- (y-mus_work)/sqrt(phi_work)
			    u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
			    u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
			Kpsim <- matrix(0,l+q,l+q)
	        Kpsim[(l+1):(l+q),(l+1):(l+q)] <- objeto$lambda.phi*objeto$M
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
				U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1)
				U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
			    D2 <- D*matrix(dgs/phi,length(phi),p)
			    D3 <- D*matrix(v2(z)/phi,length(phi),p)			
				beta_new <- beta + solve(t(D)%*%D2)%*%t(D3)%*%(y-mus)
			    W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			    psi_new <- psi +  solve(t(W_bar)%*%W_bar2 + Kpsim)%*%U_psi
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
			if(family=="Sinh-t"){
				mus_work <- objeto$mu(objeto$theta_work[1:p])[objeto$subset]			
			    gammav_work <- objeto$theta_work[(p+1):(p+l)]			
				phi_work <- exp(objeto$W%*%gammav_work)
				z_work <- (y-mus_work)/sqrt(phi_work)
				u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
	    		u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
	  	    while(tol_EM > objeto$epsilon){
				 theta <- theta_new
		  	     beta <- theta[1:p]
	    	     mus <- objeto$mu(beta)[objeto$subset]
	    	     gammav <- theta[(p+1):(p+l)]
			     phi <- exp(objeto$W%*%gammav)
				 z <- (y-mus)/sqrt(phi)			 
				 psi <- theta[(p+1):length(theta)]
				 U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1)
	    		 D <- objeto$GradD(theta)[objeto$subset,]
	 		     D2 <- D*matrix(dgs/phi,length(phi),p)
			     D3 <- D*matrix(v2(z)/phi,length(phi),p)
				 beta_new <- beta + solve(t(D)%*%D2)%*%t(D3)%*%(y-mus)
			     W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			     psi_new <- psi + solve(t(W_bar)%*%W_bar2)%*%U_psi
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
			if(family=="Sinh-t"){
			   mus_work <- objeto$mu(objeto$theta_work[1:p])[objeto$subset]
			   h_work <- objeto$theta_work[(p+1):length(objeto$theta_work)]			
			   phi_work <- exp(objeto$N%*%h_work)			
			   z_work <- (y-mus_work)/sqrt(phi_work)
			   u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
			   u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
			Kpsim <- matrix(0,l+q,l+q)
	        Kpsim[(l+1):(l+q),(l+1):(l+q)] <- objeto$lambda.phi*objeto$M
			while(tol_EM > objeto$epsilon){
				theta <- theta_new
				beta <- theta[1:p]
			    mus <- objeto$mu(beta)[objeto$subset]
				h <- theta[(p+1):length(theta)]
			    phi <- exp(objeto$N%*%h)
				z <- (y-mus)/sqrt(phi)
				psi <- theta[(p+1):length(theta)]
				U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
				D <- objeto$GradD(theta)[objeto$subset,]
			    D2 <- D*matrix(dgs/phi,length(phi),p)
			    D3 <- D*matrix(v2(z)/phi,length(phi),p)			
				beta_new <- beta + solve(t(D)%*%D2)%*%t(D3)%*%(y-mus)
			    W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			    psi_new <- psi +  solve(t(W_bar)%*%W_bar2 + Kpsim)%*%U_psi
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
	    D <- objeto$GradD(objeto$theta_work)
		D_bar <- cbind(D,objeto$Nm)
		Kpqm <- matrix(0,p+qm,p+qm)
	    Kpqm[(p+1):(p+qm),(p+1):(p+qm)] <- objeto$lambda.mu*objeto$Mm
		if(l>0){
		  if(q>0){
			if(family=="Sinh-t"){
				mus_work <- objeto$mu(objeto$theta_work[1:p])
				hm_work <- objeto$theta_work[(p+1):(p+qm)]
			    gammav_work <- objeto$theta_work[(p+qm+1):(p+qm+l)]			
				h_work <- objeto$theta_work[(p+l+qm+1):length(objeto$theta_work)]
				phi_work <- exp(objeto$W%*%gammav_work + objeto$N%*%h_work)			
			    z_work <- (y-mus_work-objeto$Nm%*%hm_work)/sqrt(phi_work)
			    u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
			    u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
			Kpsim <- matrix(0,l+q,l+q)
	        Kpsim[(l+1):(l+q),(l+1):(l+q)] <- objeto$lambda.phi*objeto$M
			while(tol_EM > objeto$epsilon){
				theta <- theta_new
			  	beta <- theta[1:p]
			    mus <- objeto$mu(beta)[objeto$subset]
				hm <- objeto$theta_work[(p+1):(p+qm)]			
		    	gammav <- theta[(p+qm+1):(p+qm+l)]
			    h <- theta[(p+qm+l+1):length(theta)]
				phi <- exp(objeto$W%*%gammav + objeto$N%*%h)
				z <- (y-mus-objeto$Nm%*%hm)/sqrt(phi)			
				psi <- theta[(p+qm+1):length(theta)]
				U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1)
				U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
			    D2 <- D_bar*matrix(dgs/phi,length(phi),p+qm)
			    D3 <- D_bar*matrix(v2(z)/phi,length(phi),p+qm)			
				beta_new <- theta[1:(p+qm)] + solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D3)%*%(y-mus-objeto$Nm%*%hm) - rbind(matrix(0,p,1),objeto$lambda.mu*objeto$Mm%*%hm))
			    W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			    psi_new <- psi +  solve(t(W_bar)%*%W_bar2 + Kpsim)%*%U_psi
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
			if(family=="Sinh-t"){
				mus_work <- objeto$mu(objeto$theta_work[1:p])
			    hm_work <- objeto$theta_work[(p+1):(p+qm)]						
			    gammav_work <- objeto$theta_work[(p+qm+1):(p+qm+l)]			
				phi_work <- exp(objeto$W%*%gammav_work)
				z_work <- (y-mus_work-objeto$Nm%*%hm_work)/sqrt(phi_work)
				u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
	    		u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
	  	    while(tol_EM > objeto$epsilon){
				 theta <- theta_new
		  	     beta <- theta[1:p]
	    	     mus <- objeto$mu(beta)
				 hm <- theta[(p+1):(p+qm)]						
	    	     gammav <- theta[(p+qm+1):(p+qm+l)]
			     phi <- exp(objeto$W%*%gammav)
				 z <- (y-mus-objeto$Nm%*%hm)/sqrt(phi)			 
				 psi <- theta[(p+qm+1):length(theta)]
				 U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1)
	 		     D2 <- D_bar*matrix(dgs/phi,length(phi),p+qm)
			     D3 <- D_bar*matrix(v2(z)/phi,length(phi),p+qm)			
				 beta_new <- theta[1:(p+qm)] + solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D3)%*%(y-mus-objeto$Nm%*%hm) - rbind(matrix(0,p,1),objeto$lambda.mu*objeto$Mm%*%hm))
			     W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			     psi_new <- psi + solve(t(W_bar)%*%W_bar2)%*%U_psi
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
			if(family=="Sinh-t"){
			   mus_work <- objeto$mu(objeto$theta_work[1:p])
			   hm_work <- objeto$theta_work[(p+1):(p+qm)]						
			   h_work <- objeto$theta_work[(p+qm+1):length(objeto$theta_work)]			
			   phi_work <- exp(objeto$N%*%h_work)			
			   z_work <- (y-mus_work-objeto$Nm%*%hm_work)/sqrt(phi_work)
			   u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
			   u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
			Kpsim <- matrix(0,l+q,l+q)
	        Kpsim[(l+1):(l+q),(l+1):(l+q)] <- objeto$lambda.phi*objeto$M
			while(tol_EM > objeto$epsilon){
				theta <- theta_new
				beta <- theta[1:p]
			    mus <- objeto$mu(beta)
	            hm <- theta[(p+1):(p+qm)]
				h <- theta[(p+qm+1):length(theta)]
			    phi <- exp(objeto$N%*%h)
				z <- (y-mus-objeto$Nm%*%hm)/sqrt(phi)
				psi <- theta[(p+qm+1):length(theta)]
				U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			    D2 <- D_bar*matrix(dgs/phi,length(phi),p+qm)
			    D3 <- D_bar*matrix(v2(z)/phi,length(phi),p+qm)			
			    beta_new <- theta[1:(p+qm)] + solve(t(D_bar)%*%D2 + Kpqm)%*%(t(D3)%*%(y-mus-objeto$Nm%*%hm) - rbind(matrix(0,p,1),objeto$lambda.mu*objeto$Mm%*%hm))
			    W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			    psi_new <- psi +  solve(t(W_bar)%*%W_bar2 + Kpsim)%*%U_psi
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
	if(qm>0){
		D_bar <- objeto$Nm
		if(l>0){
		  if(q>0){
			if(family=="Sinh-t"){
				hm_work <- objeto$theta_work[(p+1):(p+qm)]
			    gammav_work <- objeto$theta_work[(p+qm+1):(p+qm+l)]			
				h_work <- objeto$theta_work[(p+l+qm+1):length(objeto$theta_work)]
				phi_work <- exp(objeto$W%*%gammav_work + objeto$N%*%h_work)			
			    z_work <- (y-objeto$Nm%*%hm_work)/sqrt(phi_work)
			    u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
			    u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
			Kpsim <- matrix(0,l+q,l+q)
	        Kpsim[(l+1):(l+q),(l+1):(l+q)] <- objeto$lambda.phi*objeto$M
			while(tol_EM > objeto$epsilon){
				theta <- theta_new
				hm <- objeto$theta_work[(p+1):(p+qm)]			
		    	gammav <- theta[(p+qm+1):(p+qm+l)]
			    h <- theta[(p+qm+l+1):length(theta)]
				phi <- exp(objeto$W%*%gammav + objeto$N%*%h)
				z <- (y-objeto$Nm%*%hm)/sqrt(phi)			
				psi <- theta[(p+qm+1):length(theta)]
				U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1)
				U_psi[(l+1):length(psi)] <- U_psi[(l+1):length(psi)] - objeto$lambda.phi*objeto$M%*%h
			    D2 <- D_bar*matrix(dgs/phi,length(phi),p+qm)
			    D3 <- D_bar*matrix(v2(z)/phi,length(phi),p+qm)			
				beta_new <- theta[1:(p+qm)] + solve(t(D_bar)%*%D2 + objeto$lambda.mu*objeto$Mm)%*%(t(D3)%*%(y-objeto$Nm%*%hm) - objeto$lambda.mu*objeto$Mm%*%hm)
			    W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			    psi_new <- psi +  solve(t(W_bar)%*%W_bar2 + Kpsim)%*%U_psi
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
			if(family=="Sinh-t"){
			    hm_work <- objeto$theta_work[(p+1):(p+qm)]						
			    gammav_work <- objeto$theta_work[(p+qm+1):(p+qm+l)]			
				phi_work <- exp(objeto$W%*%gammav_work)
				z_work <- (y-objeto$Nm%*%hm_work)/sqrt(phi_work)
				u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
	    		u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
	  	    while(tol_EM > objeto$epsilon){
				 theta <- theta_new
		  	     beta <- theta[1:p]
				 hm <- theta[(p+1):(p+qm)]						
	    	     gammav <- theta[(p+qm+1):(p+qm+l)]
			     phi <- exp(objeto$W%*%gammav)
				 z <- (y-objeto$Nm%*%hm)/sqrt(phi)			 
				 psi <- theta[(p+qm+1):length(theta)]
				 U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1)
	 		     D2 <- D_bar*matrix(dgs/phi,length(phi),p+qm)
			     D3 <- D_bar*matrix(v2(z)/phi,length(phi),p+qm)			
				 beta_new <- theta[1:(p+qm)] + solve(t(D_bar)%*%D2 + objeto$lambda.mu*objeto$Mm)%*%(t(D3)%*%(y-objeto$Nm%*%hm) - objeto$lambda.mu*objeto$Mm%*%hm)
			     W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			     psi_new <- psi + solve(t(W_bar)%*%W_bar2)%*%U_psi
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
			if(family=="Sinh-t"){
			   hm_work <- objeto$theta_work[(p+1):(p+qm)]						
			   h_work <- objeto$theta_work[(p+qm+1):length(objeto$theta_work)]			
			   phi_work <- exp(objeto$N%*%h_work)			
			   z_work <- (y-objeto$Nm%*%hm_work)/sqrt(phi_work)
			   u_work <- (objeto$xi[2] + 1)/(objeto$xi[2] + 4*sinh(z_work)*sinh(z_work)/(alpha^2))
			}else{
			   u_work <- rep(1,length(y))
			}
			v2 <- function(z){
		          4*sinh(z)*cosh(z)*u_work/(alpha^2*z) - tanh(z)/z
		    }
			dgs <- dg2(alpha/sqrt(u_work))
			fgs <- fg2(alpha/sqrt(u_work))
			Kpsim <- matrix(0,l+q,l+q)
	        Kpsim[(l+1):(l+q),(l+1):(l+q)] <- objeto$lambda.phi*objeto$M
			while(tol_EM > objeto$epsilon){
				theta <- theta_new
	            hm <- theta[(p+1):(p+qm)]
				h <- theta[(p+qm+1):length(theta)]
			    phi <- exp(objeto$N%*%h)
				z <- (y-objeto$Nm%*%hm)/sqrt(phi)
				psi <- theta[(p+qm+1):length(theta)]
				U_psi <- (1/2)*t(W_bar)%*%(v2(z)*z^2-1) - objeto$lambda.phi*objeto$M%*%h
			    D2 <- D_bar*matrix(dgs/phi,length(phi),p+qm)
			    D3 <- D_bar*matrix(v2(z)/phi,length(phi),p+qm)			
			    beta_new <- theta[1:(p+qm)] + solve(t(D_bar)%*%D2 + objeto$lambda.mu*objeto$Mm)%*%(t(D3)%*%(y-objeto$Nm%*%hm) - objeto$lambda.mu*objeto$Mm%*%hm)
			    W_bar2 <- W_bar*matrix((fgs-1)/4,length(phi),ncol(W_bar))
			    psi_new <- psi +  solve(t(W_bar)%*%W_bar2 + Kpsim)%*%U_psi
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
	theta_new
}
