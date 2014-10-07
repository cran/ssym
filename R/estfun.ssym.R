estfun.ssym <-
function(object){
list(mu=object$score[1:(object$p+object$qm)], phi=object$score[(object$p+object$qm+1):(length(object$score))])}
