envelope <-
function(object, reps, conf){

if(object$censored==TRUE) stop("Models fitted under censored data are not supported!!",call.=FALSE)
if(missingArg(reps)) reps <- 25 else reps <- max(25,floor(reps))
if(missingArg(conf)) conf <- 0.95

if(reps!=floor(reps) | reps<=0) stop("reps must be a positive integer!!",call.=FALSE)
if(conf<=0  | conf>=1) stop("the conf argument must be within the interval (0, 1)!!",call.=FALSE)

objeto <- object
objeto$maxiter <- 1000
objeto$epsilon <- 0.0000001
theta0 <- c(object$theta.mu,object$theta.phi)
family <- object$family
xi <- object$xi
ancho <- 2:reps
dtrm <- matrix(0,length(object$y),reps)
dtrp <- matrix(0,length(object$y),reps)
i <- 1
bar <- txtProgressBar(min=0, max=reps, initial=0, width=50, char="+", style=3)
while(i <= reps){
	ysim <- object$mu.fitted  + sqrt(object$phi.fitted)*rvgs(length(object$y),family,xi)
	objeto$y <- ysim
	if(family=="Sinh-t" | family=="Sinh-normal"){
		vP <- try(itpE2(theta0,objeto), silent=TRUE)}
	else{if(family=="Powerexp" && xi < 0){vP <- itpE3(theta0,objeto)}
         else{vP <- try(itpE(theta0,objeto), silent=TRUE)}
    }
	if(is.vector(vP)){
	  thetam <- vP[1:(object$p+sum(object$qm))]
	  thetap <- vP[(object$p+sum(object$qm)+1):(object$p+sum(object$qm)+object$l+sum(object$q))]
	  if(object$orig=="linear"){
	      mu_es <- object$pspm%*%thetam
	  }else{mu_es <- object$mu(thetam)}
      phi_es <- exp(object$pspp%*%thetap)
	  z_es <- (ysim - mu_es)/sqrt(phi_es)
	  dtrm[,i] <- sort(sqrt(object$deviance.mu.f(z_es))*ifelse(z_es>=0,1,-1))
	  dtrp[,i] <- sort(sqrt(object$deviance.phi.f(z_es))*ifelse(z_es>=0,1,-1))
	  i <- i + 1
      Sys.sleep(0.5);
      setTxtProgressBar(bar,i)
	}
}
lim <- apply(dtrm, 1, quantile, prob=(1-conf)/2)
lsm <- apply(dtrm, 1, quantile, prob=(1-(1-conf)/2))
mm <- apply(dtrm, 1, mean)
lip <- apply(dtrp, 1, quantile, prob=(1-conf)/2)
lsp <- apply(dtrp, 1, quantile, prob=(1-(1-conf)/2))
mp <- apply(dtrp, 1, mean)
close(bar)
cat("\n")
par(mfrow=c(1,2))

rmus <- residuals(object)$mu
faixa <- range(rmus,lim,lsm)
par(pty="s")
qqnorm(rmus,xlab="Quantile N(0,1)",ylab="Deviance-type residuals", ylim=faixa, type="p", main="", cex=0.3, lwd=3)
par(new=TRUE)
qqnorm(lim,axes=FALSE,xlab="",ylab="",type="l",ylim=faixa,lty=1, main="")
par(new=TRUE)
qqnorm(lsm,axes=FALSE,xlab="",ylab="", type="l",ylim=faixa,lty=1, main="")
par(new=TRUE)
qqnorm(mm,axes=FALSE,xlab="", ylab="", type="l",ylim=faixa,lty=2, main="Median submodel")

rphis <- residuals(object)$phi
faixa <- range(rphis,lip,lsp)
par(pty="s")
qqnorm(rphis,xlab="Quantile N(0,1)",ylab="Deviance-type residuals", ylim=faixa, type="p", main="", cex=0.3, lwd=3)
par(new=TRUE)
qqnorm(lip,axes=FALSE,xlab="",ylab="",type="l",ylim=faixa,lty=1, main="")
par(new=TRUE)
qqnorm(lsp,axes=FALSE,xlab="",ylab="", type="l",ylim=faixa,lty=1, main="")
par(new=TRUE)
qqnorm(mp,axes=FALSE,xlab="", ylab="", type="l",ylim=faixa,lty=2, main="Skewness/Dispersion submodel")

#list(lim=lim,lsm=lsm,lip=lip,lsp=lsp,mm=mm,mp=mp)
}
