\name{processX}
\alias{processX}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Process X matrix
}
\description{
	Rescales and transforms the X matrix according to the desired parameters, and sets all the options required by the test.
}
\usage{
processX(X, family = gaussian, alpha, intercept = T, group.sizes = rep(1, ncol(X)), 
A = ncol(X), LAD = F, composite = T, M = min(10000, max(1000, 1e+10/nrow(X)/ncol(X))))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
    input matrix, of dimension n x p; each row is an observation vector.
}
  \item{family}{
	response type (see above). Default is \code{gaussian}.
}
  \item{alpha}{
	alpha for quantile rescaling; if alpha=0, then no rescaling.
}
  \item{intercept}{
	should intercept(s) be fitted (default=TRUE) or set to zero (FALSE).
}
  \item{group.sizes}{
	the vector of group sizes for affine group lasso. The number of elements is L and sum(group.sizes) should be equal to P. If L==P, then the lasso test is employed, otherwise group lasso. Default is no groups, so \code{rep(1,ncol(X))}.
}
  \item{A}{
	if A is a matrix it tests A beta = c. If A is a vector, then it gives the indexes of the parameters to be tested. Used if family=\code{gaussian}. Default is to test beta=0, so A=ncol(X).
}
  \item{LAD}{
	set TRUE if LAD lasso test. Default is FALSE
}
  \item{composite}{
	set TRUE if composite test (O & +). Default is TRUE
}
  \item{M}{
	number of Monte Carlo Simulations to estimate the distribution \eqn{\Lambda}.
}
}

\value{
	an object containing all the variables corresponding to the rescaling and test options.
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
		alpha=0.05
		X=matrix(rnorm(N*P),N,P)
		outrescale=processX(X,gaussian,alpha)
		M=100 #Leave the default or select higher value for better level.
	
	#when H0 is not rejected
		beta_scal=0
		beta=c(rep(beta_scal, s), rep(0, P-s)) 
		y=X\%*\%beta+rnorm(N)
		out=affinelassotest(y,X,gaussian,alpha,M=M,outrescale=outrescale)
		print(out$rejectH0)
	
	#when H0 is rejected
		beta_scal=10
		beta=c(rep(beta_scal, s), rep(0, P-s)) 
		y=X\%*\%beta+rnorm(N)
		out=affinelassotest(y,X,gaussian,alpha,M=M,outrescale=outrescale)
		print(out$rejectH0)
}

