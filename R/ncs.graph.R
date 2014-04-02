ncs.graph <-
function(x,g,gam,cont){

	gam <- rbind(0,gam,0)
	t <- seq(min(x),max(x),length=cont)
	nodo <- 1
	spl <- matrix(0,cont,1)
	
	for(i in 1:cont){
		if(t[i] > x[nodo + 1]){
		  nodo <- nodo + 1
		}
		h <- x[nodo+1] - x[nodo]
		spl[i] <- ((t[i]-x[nodo])*g[nodo+1] + (x[nodo+1]-t[i])*g[nodo])/h -
		          (t[i]-x[nodo])*(x[nodo+1]-t[i])*((1 + (t[i]-x[nodo])/h)*gam[nodo+1]
				  + (1 + (x[nodo+1]-t[i])/h)*gam[nodo])/6
	}
	cbind(t,spl)
}
