lambdaqut <-
function(y,X,family=gaussian,alpha.level=0.05,M=1000,qut.standardize=TRUE,intercept=TRUE,no.penalty=NULL,offset=NULL,bootstrap=TRUE,beta0=NA,method='lasso',fixbeta0=FALSE){
#FUNCTION TO OBTAIN THE QUANTILE UNIVERSAL THRESHOLD
	
	#Hidden function for vectorized bootstrapping
	qutbootstrap <- function(curz){
		ind=sample(1:n,n,replace=TRUE)
		curX=X[ind,]
		curmuhat=curz[1:n]
		curz=curz[-(1:n)]
		if(!intercept&sum(O!=0)!=0&is.null(no.penalty)){
			curmuhat=muhatZ[ind,1]
		}
		else if(!is.null(no.penalty)|sum(O!=0)!=0){
			if(is.na(beta0)){  #if beta was estimated initially
				curbeta0=glm.fit(y=curz,x=as.matrix(A[ind,]),intercept=FALSE,family=family,offset=O[ind])$coef
				curmuhat=family$linkinv(as.matrix(A[ind,])%*%curbeta0+O[ind]) 
			}
			else curmuhat=family$linkinv(as.matrix(A[ind,])%*%beta0+O[ind]) #if beta0 was specified initially
		}
		
		#zero thresholding function
		curbp=bpfunc(X=curX,z=curz,muhatZ=curmuhat,method=method)
		return(as.vector(curbp))
	}
	
	#Hidden function for vectorized computation of the intercept for each of the Monte Carlo simulations of the NULL model
	glm0=function(y,x,family,offset=offset){
		if(!fixbeta0){ #estimate beta0 for each monte carlo
			out=glm.fit(y=y,x=x,intercept=FALSE,family=family,offset=offset)
			outbeta0=out$coefficients
		}
		else outbeta0=beta0 #if beta0 was specified initially and it is not desired to estimate for each monte carlo
		return(outbeta0)
	}
	
	#Hidden function for X'z or X'z/||z||
	bpfunc=function(X,z,muhatZ,method){

		if(method=='sqrtlasso'){
			#flare standardization
			xm = matrix(rep(colMeans(X), n), nrow = n, ncol = p, byrow = TRUE)
			x1 = X - xm
			sdxinv = 1/sqrt(colSums(x1^2)/(n - 1))
			xx = x1 * matrix(rep(sdxinv, n), nrow = n, ncol = p, byrow = TRUE)
			
			Eps=z-as.matrix(muhatZ)
			xxx=colSums(Eps^2)   
			unitEps = t(t(Eps) / sqrt(xxx))
			bp= 1/sqrt(n) * abs(t(xx) %*% unitEps)
		}
		else bp=abs(t(X)%*%(z-muhatZ))

		return(bp)
	}
	
	family=family()
	
	#Check for warnings in the case of Square Root Lasso
	if(method=='sqrtlasso'&family$family!='gaussian'){
		stop('Square root lasso is just for Gaussian family')
	}
	
	n=nrow(X)
	p=ncol(X)
	X0=X
	if(is.null(offset)) O=rep(0,n)
	else O=offset

	#Check for warnings
	if(alpha.level>1|alpha.level<0){
		warning("alpha.level is not in [0,1] interval; set to 0.05")
		alpha.level=0.05
	}

	if(M<=0){
		warning("M is <=0; set to 1000")
        M=1000
	}
	if (is.null(p) | (p <= 1)) stop("X should be a matrix with 2 or more columns")
	if (n!=length(y)) stop("Number of observations in y not equal to number of rows in X")
	if(length(O)!=n) stop("length of offset is different to the number of observations")
	if(family$family=='poisson'&(sum(y==0)==n)) stop("All your Poisson counts are zero")
	if(family$family=='binomial'&((sum(y==0)==n)|(sum(y==1)==n))) stop("All your Binomial measures are zero or one")
	
	#initialize A matrix
	if(!intercept&is.null(no.penalty)){
		muhat=family$linkinv(O)
		muhatfull=muhat
	}
	else{
		A=c()
		if(!is.null(no.penalty)){ 
			A=as.matrix(X[,no.penalty]) #if there are more unpenalized coefficients
			X=X[,-no.penalty]
		}
		if(intercept) A=cbind(rep(1,n),A)  #if there is an intercept (column of ones)
		
		#glm(y~A) required for obtaining lambda.max
		beta0full=glm.fit(y=y,x=as.matrix(A),intercept=FALSE,family=family,offset=offset)$coef
		
		if(!is.numeric(beta0)){
			#Estimate beta0 as glm(y~A)
			beta0=beta0full
		}
		else if(length(beta0)!=(length(no.penalty)+intercept)) stop("length of beta0 is different from the number of unpenalized covariates or the intercept has not been included")
		muhat=family$linkinv(as.matrix(A)%*%beta0+O)
		muhatfull=family$linkinv(as.matrix(A)%*%beta0full+O)
	}
		
	#Set default alpha.level
	if(alpha.level=='default')	alpha.level=1/(sqrt(pi*log(p)))

	#Check if A=1 and there is an explicit characterization of D
	if(!(intercept&is.null(no.penalty)&sum(O==0)==n)) warning("Explicit characterization of D is not defined")

	#Do valid Monte Carlo Simulation of the null model
	znotinD=0

	for(nMC in 1:2){
		if(family$family=='gaussian') z=matrix(rnorm(n*M,mean=as.vector(muhat),sd=1),n,M)
		else if(family$family=='poisson'){
			z=matrix(rpois(n*M,lambda=as.vector(muhat)),n,M)
			znotinD=sum(apply(z,2,sum)==0)
			z=z[,apply(z,2,sum)!=0]
		}	
		else if(family$family=='binomial'){
			z=matrix(rbinom(n*M,size=1,prob=as.vector(muhat)),n,M)
			znotinD=sum((apply(z,2,sum)==0)|(apply(z,2,sum)==n))
			z=z[,(apply(z,2,sum)!=0)&(apply(z,2,sum)!=n)]
		}
		else stop("Not available family")
		
		#Check if it was a valid simulation
		if(znotinD<=(M-2)) break
		else{ #NOT VALID MC simulation
			#reestimate the intercept without interation
			warning('Intercept will be estimated without iteration since there were no valid simulations under null hypothesis with current intercept')
			beta0=glm.fit(y=y,x=as.matrix(A),intercept=FALSE,family=family,offset=offset)$coef
			muhat=family$linkinv(as.matrix(A)%*%beta0+O)
		}
	}
	if(znotinD>(M-2)) stop("Can't generate valid simulations under the null hypothesis, try increasing M or changing the intercept")
	
	#obtain estimate of beta0 for each simulation
	if(!intercept&is.null(no.penalty)){
		muhatZ=family$linkinv(O)
		muhatZ=matrix(rep(muhatZ,ncol(z)),nrow=n)
	}
	else{
		beta0Z=apply(z,2,glm0,x=as.matrix(A),family=family,offset=offset)
		muhatZ=family$linkinv(as.matrix(A)%*%beta0Z+O)
	}

	divX=rep(1,ncol(X))
	
	if(!bootstrap){
		#No Bootstrap Matrix X
		
		#zero thresholding function

		bp=bpfunc(X=X,z=z,muhatZ=muhatZ,method=method)
		scale.factor=rep(1,p)
		
		#quantile-based standardization
		if(qut.standardize){
			divX=apply(abs(bp),1,quantile,prob=(1-alpha.level))
			divX[which(divX==0)]=1
			
			#Alignment/scaling
			X=t(t(X)/divX)
			
			#zero thresholding function for the standardized X matrix
			bp=bpfunc(X=X,z=z,muhatZ=muhatZ,method=method)
		}
	}
	else{
		#Bootstrap Matrix X

		#quantile-based standardization
		if(qut.standardize&method!='sqrtlasso'){
		
			#zero thresholding function
					
			bp=bpfunc(X=X,z=z,muhatZ=muhatZ,method=method)
			scale.factor=rep(1,p)
			divX=apply(abs(bp),1,quantile,prob=(1-alpha.level))
			divX[which(divX==0)]=1
			
			#Alignment/scaling
			X=t(t(X)/divX)
		}
		#Bootstrapping
		bp=apply(rbind(muhatZ,z),2,qutbootstrap)
	}
	
	#Obtain distribution of Lambda=max(X'z)
	resultsMC=apply(bp,2,max,na.rm=TRUE)
	resultsMC=c(resultsMC,rep(Inf,znotinD))
	#obtain lambda.qut for the given quantile
	lambda=quantile(resultsMC, prob=(1-alpha.level))
	
	#lambda max
	#bp=abs(t(X)%*%(y-muhat))
	bp=bpfunc(X=X,z=y,muhatZ=muhatfull,method=method)
	lambdamedian=median(bp)
	lambdamax=max(bp)
	
	
	if(!is.null(no.penalty)){
		X0[,-no.penalty]=X
		scale.factor[-no.penalty]=divX
		X=X0
	}
	else scale.factor=divX
		
	#OUTPUT
	out=NULL
	out$scale.factor=scale.factor
	out$lambda.max=lambdamax
	out$lambda.median=lambdamedian
	out$lambda=lambda
	out$beta0=beta0
	out$Xnew=X
	return(out)
	
}
