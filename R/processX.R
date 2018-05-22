processX=function(X,family=gaussian,alpha,intercept=T,group.sizes=rep(1,ncol(X)),A=ncol(X), LAD=F, composite=T, M=min(1.e4, max(1000,1.e10/nrow(X)/ncol(X)))){
  ## linear model y=intercept+X beta+epsilon
  ##    s.t. A beta = b, if A is a matrix 
  ##    s.t. beta1=b, where beta=(beta0, beta1) and length(beta1)=A, if A is an integer
  ## Note that b does not play a role for quantile rescaling (I-P_{K_A}) XA'(AA')^{-1}
  ## alpha for quantile rescaling; if alpha=0, then no rescaling

  L=length(group.sizes)
  N=nrow(X); P=ncol(X)
  X0=NA;
  
  c0=(family()$family!="gaussian"& is.matrix(A))||(family()$family!="gaussian"& LAD)
  c1=(LAD & is.matrix(A))
  if(LAD || family()$family!="gaussian") {X0=X[,0:(P-A), drop=F]; if(intercept) X0=cbind(1, X0)}
  
  if(!is.matrix(A)){A=cbind(matrix(0,nrow=A,ncol=ncol(X)-A) ,diag(A))}
  
  nrowA=nrow(A)
  qrA=qr(A)
  c2=(qrA$rank != nrowA)
  c3=(sum(group.sizes)!=nrowA)
  if(c0 || c1 || c2 || c3) {
    if(c0) stop("Affine lasso or LAD is currently only for Gaussian, A should not be a matrix")
    if(c1) stop("Affine LAD not available. LAD lasso and LAD group lasso are.")
    if(c2) stop("The matrix A in H0: A beta = b is not full row rank!")
    if(c3) stop("The group.sizes do not match the nbr of constraint in H0.")
    rankXkerA=NA
    lambdas.alpha=NA
    XkerAProj=NA
    XkerAProj.mat=NA
    X0=NA
  }
  else{
	  X=cbind(intercept, X)
	  A=cbind(!intercept, A)
	  B=X %*% t(A) %*% solve(A%*%t(A))
	  Av=svd(A, nv=ncol(A))$v
	  if(ncol(A)>nrowA) { ## A has a kernel > {0}
		K_A=Av[,(nrowA+1):ncol(A), drop=F]
		XK_A=X %*% K_A
		qrXK_A=qr(XK_A)
		rankXkerA=qrXK_A$rank
		if(rankXkerA<N){
		  XkerAProj.mat=XK_A %*% solve(t(XK_A) %*% XK_A) %*% t(XK_A)
		  XkerAProj = function(y, XkerAProj.mat){
			XkerAProj.mat %*% y
 	      }
		  X = B - apply(B,2,XkerAProj, XkerAProj.mat=XkerAProj.mat)
		}
	  }
	}
    
  ##############
  ## Rescaling
  if(rankXkerA<N){
    if(alpha>0){
      if(max(group.sizes)==1){
        ##rescaling by sd
        lambdas.alpha=apply(X,2,sd)
        X=t(t(X)/lambdas.alpha)
		sf=lambdas.alpha
      }
      else {
        # Blockwise quantile rescaling
        lambdas.alpha=rep(NA,L)
        i=0
        Xr=c()
        sf=c() #scale factor
        z=matrix(rnorm(N*M), N, M)
        if(LAD && (ncol(X0)>0)) {outlm=lm(z~X0-1); z=z-predict(outlm)}
        for(l in 1:L){
          Xl=X[,i+1:group.sizes[l]]
          if(!LAD) u=t(Xl) %*% z
          else {
            u=t(Xl) %*% sign(z)
          }
          lambdas=sqrt(apply(u^2,2,sum))
          lambdas.alpha[l]=quantile(lambdas, 1-alpha)
          Xr=cbind(Xr, Xl/lambdas.alpha[l])
          sf=c(sf,rep(lambdas.alpha[l],group.sizes[l]))
          i=i+group.sizes[l]
        }
        
        if(composite){
          lambdas.alpha.composite=apply(X,2,sd)
          Xr2=t(t(X)/lambdas.alpha.composite)
          lambdas=apply(abs(t(Xr2) %*% z), 2, max)
          lambda0=quantile(lambdas, 1-alpha)
          Xr2=Xr2/lambda0
          lambdas.alpha.composite=lambdas.alpha.composite*lambda0
          sf2=lambdas.alpha.composite
        }
        X=Xr
        if(composite){
        	X=cbind(X,Xr2)
        	lambdas.alpha=c(lambdas.alpha, lambdas.alpha.composite)
        	sf=c(sf,sf2)
        }
      }
    }
    else {
      lambdas.alpha=rep(1,L)
      sf=rep(1,P)
    }
  }
  else {
    stop("X ker(A) is a full row rank matrix: your null hypothesis cannot be tested!")
  }
 
  out=NULL
  out$rankXkerA=rankXkerA
  out$lambdas.alpha=lambdas.alpha
  out$Xrescaled=X
  out$sf=sf
  out$XkerAProj=XkerAProj
  out$XkerAProj.mat=XkerAProj.mat
  out$group.sizes=group.sizes
  out$LAD=LAD
  out$X0=X0
  out$composite=composite
  return(out)
}
