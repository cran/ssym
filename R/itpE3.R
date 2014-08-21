itpE3 <-
function(vP,objeto){
		 objeto$theta_work <- vP
		 theta_new <- itpM3(objeto)
	list(theta=theta_new,criterion=theta_new,iter=1)
}
