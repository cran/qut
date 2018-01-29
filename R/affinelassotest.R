affinelassotest=function(y,Xdata,family=gaussian,alpha,cc=NA,lambda.alpha=NA,outrescale=NA,intercept=T,group.sizes=rep(1,ncol(X)),A=ncol(X), LAD=F, composite=T, M=round(min(1.e4, max(1000,1.e9/nrow(X)/ncol(X))))){
  ##
  ## Test H0: A beta = cc
  ##
  ## o y=response vector of length N
  ## o X=matrix of covariates of dimension NxP
  ## o alpha=level of the test
  ## o intercept: if TRUE, then an intercept is included (and not tested)
  ## o group.sizes is the vector of group sizes for affine group lasso. The number of elements is L and sum(group.sizes) should be equal to P. If L==P, then the lasso test is employed, otherwise group lasso
  ## o if A is a matrix and it is affine lasso. If A is a vector, then it gives the indexes of the parameters to be tested.
  ## This is to test H=0: beta[subset given in A]=0 directly by using directly the simple formulae of the zero thresholding function of lasso

  X=Xdata

  N=length(y); P=ncol(X);
  if(nrow(X) != N) {
    print("The number of rows of X is different from the length of y")
  }
  else{
    if(sum(!is.na(outrescale))==0){
      outrescale=processX(X,family=family,alpha,intercept,group.sizes,A,LAD,composite,M)
    }
    Xscaled=outrescale$Xrescaled
    XkerAProj=outrescale$XkerAProj
    XkerAProj.mat=outrescale$XkerAProj.mat
    rankXkerA=outrescale$rankXkerA
    group.sizes=outrescale$group.sizes
    LAD=outrescale$LAD
    composite=outrescale$composite
    
    if(sum(!is.na(lambda.alpha))==0){
      if(family()$family=="gaussian") z=matrix(rnorm(N*M), N, M)
      if(family()$family=="poisson") z=matrix(rpois(N*M, lambda=1), N, M)
      if(family()$family=="binomial") z=matrix(rbinom(N*M, 1, prob=0.5), N, M)
      
      lambdas=apply(z,2,ztf, Xdata=X,family=family,intercept=intercept, group.sizes=group.sizes, LAD=LAD, outrescale=outrescale,composite=composite)
      lambda.alpha=quantile(lambdas, 1-alpha)
    }
    else{lambdas=NA}
    lambda.data=ztf(y,Xdata=X,family=family,intercept=intercept, group.sizes=group.sizes, LAD=LAD, outrescale=outrescale,composite=composite)
    rejectH0=(lambda.data>lambda.alpha)
  } 
       
  out=NULL
  out$lambda.alpha=lambda.alpha
  out$lambda.data=lambda.data
  out$rejectH0=rejectH0
  out$lambdas=lambdas
  out$outrescale=outrescale
  return(out)
}
#out0lasso=affinelassotest(yy,X,family,alpha, group.sizes=rep(1,ncol(X)), A=A)