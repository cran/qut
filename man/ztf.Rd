\name{ztf}
\alias{ztf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zero thresholding function
}
\description{
	Obtains the value of the minimum regularization parameter that sets all coefficients to zero for different types of thresholding tests.
}
\usage{
ztf(y,Xdata, family=gaussian, A=ncol(Xdata), cc=NA, intercept=T,
group.sizes=rep(1,ncol(Xdata)), LAD=F, outrescale=NA, composite=T, alpha=0,
M=min(1.e4, max(1000,1.e10/nrow(Xdata)/ncol(Xdata))))
}

\arguments{
  \item{y}{
	response variable. Quantitative for family=\code{gaussian}, or family=\code{poisson} (non-negative counts). For family=\code{binomial} should be a factor with two levels.
}
  \item{Xdata}{
     input matrix, of dimension n x p; each row is an observation vector.
}
  \item{family}{
	response type (see above). Default is \code{gaussian}.
}
  \item{A}{
	if A is a matrix it tests A beta = c. If A is a vector, then it gives the indexes of the parameters to be tested. Used if family=\code{gaussian}. Default is to test beta=0, so A=ncol(X).
}
  \item{cc}{
	vector c
}
  \item{intercept}{
	should intercept(s) be fitted (default=TRUE) or set to zero (FALSE).
}
  \item{group.sizes}{
	the vector of group sizes for affine group lasso. The number of elements is L and sum(group.sizes) should be equal to P. If L==P, then the lasso test is employed, otherwise group lasso. Default is no groups, so \code{rep(1,ncol(X))}.
}
  \item{LAD}{
	set TRUE if LAD lasso test. Default is FALSE
}
  \item{outrescale}{
	object containing all variables corresponding to the rescaling and test options. If not provided, this is calculated authomatically with function \code{processX}. Default is NA.
}
  \item{composite}{
	set TRUE if composite test (O & +). Default is TRUE
}
  \item{alpha}{
	alpha for quantile rescaling; if alpha=0, then no rescaling.
}
  \item{M}{
	number of Monte Carlo Simulations to estimate the distribution \eqn{\Lambda}.
}
}
\value{
	value of the minimum regularization parameter that sets all coefficients to zero
}
\author{
	Sylvain Sardy and Jairo Diaz
}
\examples{
	# Test H0:beta=0
		P=200
		N=20
		s=1
		A=P
		X=matrix(rnorm(N*P),N,P)
		M=100 #Leave the default or select higher value for better level.
	#when H0 is not rejected
		beta_scal=0
		beta=c(rep(beta_scal, s), rep(0, P-s)) 
		y=X\%*\%beta+rnorm(N)
		zerolambda=ztf(y,X,M=M)
		print(zerolambda)
	
	#when H0 is rejected
		beta_scal=10
		beta=c(rep(beta_scal, s), rep(0, P-s)) 
		y=X\%*\%beta+rnorm(N)
		zerolambda=ztf(y,X,M=M)
		print(zerolambda)
}
