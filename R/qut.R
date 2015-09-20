qut <-
function(y,X,fit,family=gaussian,alpha.level='default',M=1000,qut.standardize=TRUE,intercept=TRUE,offset=NULL,bootstrap=TRUE,sigma='qut',estimator='unbiased',type=c('glmnet','lars'),lambda.seq=0,penalty.factor=rep(1,p),lambda.min.ratio=ifelse(n<p,0.01,0.0001),nlambda=100,lambda=NULL,...){
#FUNCTION TO FIT GLM WITH THE QUANTILE UNIVERSAL THRESHOLD

	#Hidden function to get glm fit with non-zero coefficients found by LASSO-GLM
	betaGLMestimates=function(){
		if(type=='glmnet') betatemp=beta[-1] else betatemp=beta
		betaglm=betatemp*0
		
		if(is.null(offset)) O=rep(0,n)
		if(sum(betatemp!=0)>0){
			if(intercept) betaglm0=glm(y~X0[,betatemp!=0],family=family,offset=offset)$coef
			else betaglm0=glm(y~X0[,betatemp!=0]-1,family=family,offset=offset)$coef
		}
		else{if(intercept) betaglm0=glm(y~1,family=family,offset=offset)$coef}
		if(intercept){
			betaglm[betatemp!=0]=betaglm0[-1]
			betaglm=c(betaglm0[1],betaglm)
		}
		else{
			betaglm[betatemp!=0]=betaglm0
			betaglm=c(0,betaglm)
		}
		
		names(betaglm)[1]='(Intercept)'	
		names(betaglm)[-1]=1:p
		return(betaglm)
	}

	#Initialize some basic variables
	X0=X
	n=nrow(X)
	p=ncol(X)
	f=family()
	no.penalty=which(penalty.factor==0)
	if(length(no.penalty)==0) no.penalty=NULL
	
	#Check for warnings
	type=match.arg(type)
	if(type=='lars'&f$family!='gaussian'){
		warning("type lars is just for gaussian family; set to glmnet")
		type='glmnet'
	}
	if(type=='lars'&!is.null(no.penalty)){
		warning("type lars does not allow no penalty subsets; set to glmnet")
		type='glmnet'
	}
	if(type=='lars'&!is.null(offset)){
		warning("type lars does not allow offset; set to glmnet")
		type='glmnet'
	}
	
    #Obtain lambda.qut (Quantile Universal Threshold)
	outqut=lambdaqut(y=y,X=X,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize, family=family,intercept=intercept,no.penalty=no.penalty,offset=offset,bootstrap=bootstrap)
	X=outqut$Xnew #Get standardized matrix
	lambdaqut=outqut$lambda
	
	#Estimate sigma for gaussian family (if not specified)
	if(f$family=='gaussian'&!is.numeric(sigma)){
		if(sigma=='qut') sigma=sigmaqut(y=y,X=X,intercept=intercept,estimator=estimator,alpha.level=alpha.level,M=M,qut.standardize=qut.standardize,offset=offset,penalty.factor=penalty.factor,...)
		else{
			warning("sigma rcv or cv does not take into account penalty factor")
			if(sigma=='rcv') sigma = sigmarcv(X=X,y=y,cv=FALSE,intercept=intercept)$sigmahat
			else if(sigma=='cv') sigma = sigmarcv(X=X,y=y,cv=TRUE,intercept=intercept)$sigmahat
			else sigma = 1
		}
	}
	else sigma=1
	
	#Fit model
	if(type=='glmnet'){ #GLMNET
		lambdamax=outqut$lambda.max*sum(penalty.factor)/p
		#set sequence of lambdas
		if(lambda.seq==2){  #default glmnet values
			if(missing(lambda.min.ratio)) lambda.min.ratio=ifelse(n<p,0.01,0.0001)
			if(missing(nlambda)) nlambda=100
			if(missing(lambda)) lambda=NULL
		}
		else{ #starting from lambdaqut to lambdamax
			lambda.min.ratio=lambdaqut*sigma/lambdamax*0.9999
			if(lambda.min.ratio>1) lambda.min.ratio=0.9999
			if(lambda.seq==0) nlambda=100 #100 lambdas  #DEFAULT
			else if(lambda.seq==1) nlambda=n #n lambdas
			lambda=exp(seq(log(lambdamax),log(lambdaqut*sigma),length=nlambda))/n #n lambdas in logscale
		}
		if(missing(penalty.factor)) penalty.factor=rep(1,p)
		if(!missing(no.penalty)) penalty.factor[no.penalty]=0

		if(missing(fit)) fit=glmnet(X,y,standardize=FALSE,intercept=intercept,family=f$family,penalty.factor=penalty.factor,offset=offset,nlambda=nlambda,lambda.min.ratio=lambda.min.ratio,...)		
		
		if(max(fit$lambda)==Inf) beta=rep(0,p+1)
		else beta=coef(fit, s=lambdaqut*sigma/n*sum(penalty.factor)/p,offset=offset)/c(1,outqut$scale.factor) 

		lambdaqut=lambdaqut/n

		#Get GLM fitted coefficients
		betaglm=try(betaGLMestimates(),silent=TRUE)

		
	}
	else if(type=='lars'){ #LARS
	
		#LASSO fit
		lambdamax=outqut$lambda.max*n
		if(missing(penalty.factor)) penalty.factor=rep(1,p)
		X=X/penalty.factor
		
		if(missing(fit)) fit=lars(X,y,intercept=intercept,normalize=FALSE,...)
		beta=coef(fit,s=lambdaqut*sigma,mode='lambda')
		
		if(intercept) beta0=mean(y)-mean(X%*%beta) else beta0=0
		beta=beta/penalty.factor/outqut$scale.factor
		
		#Least squares fit
		betaglmestimatesLARS=function(){
			betaglm=beta*0
			if(sum(beta!=0)>0){
				if(intercept) betaglm0=lm(y~X0[,beta!=0])$coef
				else betaglm0=lm(y~X0[,beta!=0]-1)$coef
			}
			else{if(intercept) betaglm0=lm(y~1)$coef}
			
			if(intercept){
				betaglm[beta!=0]=betaglm0[-1]
				betaglm=c(betaglm0[1],betaglm)
			}
			else{
				betaglm[beta!=0]=betaglm0
				betaglm=c(0,betaglm)
			}
			names(betaglm)[1]='(Intercept)'	
			names(betaglm)[-1]=1:p
			return(betaglm)
		}
		
		#Get GLM fitted coefficients
		betaglm=try(betaGLMestimates(),silent=TRUE)
		
		beta=c(beta0,beta)

		names(beta)[1]='(Intercept)'	
		names(beta)[-1]=1:p	
		
	}

	#OUTPUT
	out=NULL
	out$lambda=lambdaqut
	out$lambda.max=lambdamax/n
	out$scale.factor=outqut$scale.factor
	out$beta=beta
	if(class(betaglm)[1]=="try-error") warning("No valid GLM fits were possible to find")
	out$betaglm=betaglm
	out$fit=fit
	out$family=family
	if(f$family=='gaussian') out$sigma=sigma
	else out$sigma=NA
	class(out)='qut'
	return(out)
}


