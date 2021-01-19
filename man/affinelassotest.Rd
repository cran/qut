\name{affinelassotest}
\alias{affinelassotest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Affine lasso test
}
\description{
	Perform thresholding tests for null hypothesis of the form H0: A beta = c.
}
\usage{
affinelassotest(y, Xdata, family = gaussian, alpha, beta0=NA, cc = NA, lambdas = NA, 
outrescale = NA, intercept = TRUE, group.sizes = rep(1, ncol(X)), A = ncol(X), 
LAD = FALSE, composite = TRUE, M = round(min(10000, max(1000, 1e+09/nrow(X)/ncol(X)))))
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
  \item{alpha}{
	desired level of the test.
}
  \item{beta0}{
	if known, value of the nuissance parameter. Otherwise it takes beta0=0, and the statistic is asymptotically a pivot for non-Gaussian.
}
  \item{cc}{
	vector c
}
  \item{lambdas}{
	if not provided, the code performs M Monte Carlo simulation to obtain the empirical distribution \eqn{\Lambda} and the corresponding value in the \eqn{\alpha}-quantile for testing. Otherwise, vector with all the Monte Carlo values of lambda. Default is NA.
}
  \item{outrescale}{
	object containing all variables corresponding to the rescaling and test options. If not provided, this is calculated authomatically with function \code{processX}. Default is NA.
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
  \item{lambda.alpha}{value of \eqn{\lambda} in the \eqn{\alpha}-quantile}
  \item{lambda.data}{value of \eqn{\lambda} for the current data.}
  \item{rejectH0}{result of the test. TRUE if H0 is rejected.}
  \item{lambdas}{values of \eqn{\lambda} of the Monte Carlo simulation under the null hypothesis.}
  \item{outrescale}{object containing all the rescaling variables.}
  \item{pvalue}{approximate p-value for current data.}
}
\references{
	Thresholding tests. Sylvain Sardy, Caroline Giacobino, Jairo Diaz.
}
\author{
	Sylvain Sardy and Jairo Diaz Rodriguez
}

\examples{
	# Test H0:beta=0
		P=200
		N=20
		s=1
		A=P
		alpha=0.05
		X=matrix(rnorm(N*P),N,P)
		M=100 #Leave the default or select higher value for better level.
	
	#when H0 is not rejected
		beta_scal=0
		beta=c(rep(beta_scal, s), rep(0, P-s)) 
		y=X\%*\%beta+rnorm(N)
		out=affinelassotest(y,X,gaussian,alpha,M=M)
		print(out$rejectH0)
	
	#when H0 is rejected
		beta_scal=10
		beta=c(rep(beta_scal, s), rep(0, P-s)) 
		y=X\%*\%beta+rnorm(N)
		out=affinelassotest(y,X,gaussian,alpha,M=M)
		print(out$rejectH0)
}

