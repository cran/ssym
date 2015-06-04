itpEC2 <-
function(vP,objeto){
	epsilon <- objeto$epsilon
	maxiter <- objeto$maxiter
	y <- objeto$y
	event <- objeto$event
	n <- length(y)
	pspm <- objeto$pspm
    penm <- objeto$penm
    pspp <- objeto$pspp
    penp <- objeto$penp
	p <- objeto$p
	qm <- objeto$qm
	q <- objeto$q
	l <- objeto$l
	theta_new <- vP
	v <- objeto$v
	vp <- objeto$vp	
	h <- objeto$h	
	tol <- 1
	cont <- 0
	Ltt <- matrix(0,ncol(pspm)+ncol(pspp),ncol(pspm)+ncol(pspp))
	Ut <- matrix(0,ncol(pspm)+ncol(pspp),1)
	while(tol > epsilon){
		theta <- theta_new
  		mu_es <-  pspm%*%theta[1:(p + sum(qm))]
		phi_es <- exp(pspp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)])
		z_es <- (y - mu_es)/sqrt(phi_es)
		h_es <- ifelse(event==1,h(z_es),0)
		v_es <- v(z_es)
		vp_es <- ifelse(event==1,0,vp(z_es))
		Dc   <- ifelse(event==1,h_es*(h_es-v_es*z_es),vp_es*z_es + v_es)
		Dcar <- ifelse(event==1,h_es + z_es*h_es*(h_es-v_es*z_es),vp_es*z_es^2 + 2*v_es*z_es)/2
		Dcab <- Dcar*z_es/2
		Ltt[1:ncol(pspm),1:ncol(pspm)] <- -crossprod(pspm,pspm*matrix(Dc/phi_es,n,ncol(pspm))) - penm
		Ltt[1:ncol(pspm),1:ncol(pspm)] <- -crossprod(pspm,pspm*matrix(Dc/phi_es,n,ncol(pspm))) - penm
		Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)] <- -crossprod(pspp,pspm*matrix(Dcar/sqrt(phi_es),n,ncol(pspm)))
		Ltt[1:ncol(pspm),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- t(Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),1:ncol(pspm)])
		Ltt[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp)),(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- -crossprod(pspp,pspp*matrix(Dcab,n,ncol(pspp))) - penp
		vt_es <- ifelse(event==1,h_es/z_es,v_es)
		s_es <- ifelse(event==1,h_es*z_es+1,v_es*z_es^2)		
		Ut[1:ncol(pspm)] <- crossprod(pspm,vt_es*z_es/sqrt(phi_es)) - penm%*%theta[1:(p+sum(qm))]
		Ut[(ncol(pspm)+1):(ncol(pspm)+ncol(pspp))] <- crossprod(pspp,(s_es - 1)/2) - penp%*%theta[(p+sum(qm)+1):(p+sum(qm)+sum(q)+l)]
		theta_new <- theta + crossprod(solve(-Ltt),Ut)
		theta_pd <- ifelse(abs(theta_new) <= epsilon,1,theta_new)
		tol <- max(abs(theta_new-theta)/abs(theta_pd))
		cont <- cont + 1
		if(cont > maxiter) stop("no convergence was obtained!!",call.=FALSE)
    }
	theta_new
}
