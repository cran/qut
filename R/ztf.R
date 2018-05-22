ztf=function(y,Xdata,family=gaussian,A=ncol(Xdata),cc=NA,intercept=T,group.sizes=rep(1,ncol(Xdata)),LAD=F,outrescale=NA,composite=T,alpha=0,M=min(1.e4, max(1000,1.e10/nrow(Xdata)/ncol(Xdata)))){
  ## A beta = cc
  
  if(min(y)==max(y)){
    if((family()$family=="poisson"&y[1]==0)|family()$family!="poisson"){
      stop("It is not possible to get the zero thresholding function of degenerated data")
    }
  } 
  
  if(family()$family!="gaussian"& LAD) stop("LAD is not available for family different from Gaussian")
  if(!is.na(cc)&family()$family!="gaussian") stop("Currently, it is not possible to test beta=cc when cc different from zero for family different from Gaussian")
  if(family()$family!="gaussian") warning("The A matrix should be of the form [0 I] when family is different from Gaussian")
  
  X=Xdata
  N=nrow(X)
  P=ncol(X)
  
  if(sum(!is.na(outrescale))==0){
    outrescale=processX(X,family=family,alpha,intercept,group.sizes,A, LAD, composite, M)
  }

  X0=outrescale$X0
  Xscaled=outrescale$Xrescaled
  XkerAProj=outrescale$XkerAProj
  XkerAProj.mat=outrescale$XkerAProj.mat

  L=length(group.sizes)
  if(!is.na(cc)) y=y-Xscaled %*% cc # to test whether in confidence region
  
  #center y (not necessary for gaussian)
  if(family()$family!="gaussian"){
    if((dim(X0)[2]==1)&(sum(X0)==length(y))){
      muhat=rep(mean(y),length(y))
    }
    else{
    	#warning("You might need to mean center X")
    	lmfit=glm(y~X0-1,family=family)
    	muhat=predict(lmfit, type="response")
    }
  	muhatbar=mean(muhat)
  	yCentered=y-muhat
  }
  else yCentered=y
  
  if(max(group.sizes)==1){
  # lasso
    if(LAD) {outlm=lm(y~X0-1); y=y-predict(outlm); denom=1; num=max(abs(t(Xscaled) %*% sign(y)))}
    else {

      if(family()$family=="gaussian"&(length(y)>(2*ncol(X)))){
      ## Pivot with RSS under full model
        lmfit=glm(y~cbind(intercept,X),family=family)
		denom=sqrt(N*mean(lmfit$res^2))
      }
      else {
      ## Pivot with RSS under H0, i.e., square-root
	    
		#choose the correct denominator of the test statistic
        if(family()$family=="gaussian") denom=sqrt(mean((y-XkerAProj(y, XkerAProj.mat=XkerAProj.mat))^2))  
		else if(family()$family=="poisson")	denom=sqrt(N*muhatbar)
		else if(family()$family=="binomial") denom=sqrt(N*(muhatbar*(1-muhatbar)))
      }
      num=max(abs(t(Xscaled) %*% yCentered))  
    }
	
  }
  
  
  if(max(group.sizes)!=1 & composite==F){
  ## group lasso
    if(LAD) {outlm=lm(y~X0-1); y=y-predict(outlm); denom=1}
    num=0
    i=0
    for(l in 1:L){
      Xl=Xscaled[,i+1:group.sizes[l]]
      if(!LAD) u=t(Xl) %*% yCentered
      else u=t(Xl) %*% sign(y)
      lambdas=sqrt(sum(u^2))
      num=max(num, lambdas)
      i=i+group.sizes[l]
    }
    if(!LAD){
      if(family()$family=="gaussian"&(length(y)>(2*ncol(X)))){
      ## Pivot with RSS under H1
        lmfit=lm(y~cbind(intercept,X))
		denom=sqrt(N*mean(lmfit$res^2))
      }
      else {
      ## Pivot with RSS under H0, i.e., square-root
	  
        #choose the correct denominator of the test statistic
        if(family()$family=="gaussian") denom=sqrt(mean((y-XkerAProj(y, XkerAProj.mat=XkerAProj.mat))^2))  
		else if(family()$family=="poisson")	denom=sqrt(N*muhatbar)
		else if(family()$family=="binomial") denom=sqrt(N*(muhatbar*(1-muhatbar))) 
      }
    }
  }

  #O&+ test
  if(max(group.sizes)!=1 & composite==T){
    ## composite lasso
    P=ncol(Xscaled)/2
    X1=Xscaled[,c(1:P)]
    X2=Xscaled[,-c(1:P)]
    
    if(LAD) {outlm=lm(y~X0-1); y=y-predict(outlm); denom=1}
    num=0
    i=0
    for(l in 1:L){
      Xl=X1[,i+1:group.sizes[l]]
      if(!LAD) u=t(Xl) %*% yCentered
      else u=t(Xl) %*% sign(y)
      lambdas=sqrt(sum(u^2))
      num=max(num, lambdas)
      i=i+group.sizes[l]
    }
    if(!LAD) num=max(num, max(abs(t(X2) %*% yCentered)))
    else num=max(num, max(abs(t(X2) %*% sign(y))))
    if(!LAD){
      if(family()$family=="gaussian"&(length(y)>(2*ncol(X)))){
        ## Pivot with RSS under H1
        lmfit=lm(y~cbind(intercept,X))
        denom=sqrt(N*mean(lmfit$res^2))
      }
      else {
        ## Pivot with RSS under H0, i.e., square-root

		  #choose the correct denominator of the test statistic
          if(family()$family=="gaussian") denom=sqrt(mean((y-XkerAProj(y, XkerAProj.mat=XkerAProj.mat))^2))  
		  else if(family()$family=="poisson")	denom=sqrt(N*muhatbar)
		  else if(family()$family=="binomial") denom=sqrt(N*(muhatbar*(1-muhatbar)))		
      }
    }
  }
  
  return(num/denom)
} 
