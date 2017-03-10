coef.qut <- 
function(object, mode='glm',...){
	if(mode=='lasso') return(object$beta)
	else if(class(object$betaglm)=='try-error'){
		warning('GLM did not converged, returning LASSO coefficients')
		return(object$beta)
	}
	else{
		if((!attr(object$betaglm,'converged'))){
			warning('GLM did not converged, returning LASSO coefficients')
			return(object$beta)
		}
		else return(object$betaglm)
	}
}