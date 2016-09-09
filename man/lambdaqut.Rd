\name{lambdaqut}
\alias{lambdaqut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Quantile Universal Threshold, regularization parameter for GLM-lasso
}
\description{
Computes the Quantile Universal Threshold for GLM-lasso.
}
\usage{
lambdaqut(y, X, family = gaussian, alpha.level = 0.05, M = 1000, 
qut.standardize = TRUE, intercept = TRUE, no.penalty = NULL, offset = NULL,
bootstrap=TRUE,beta0=NA,method='lasso')
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	response variable. Quantitative for family=\code{gaussian}, or family=\code{poisson} (non-negative counts). For family=\code{binomial} should be a factor with two levels.
}
  \item{X}{
     input matrix, of dimension n x p; each row is an observation vector.
}
  \item{family}{
	response type (see above). Default is \code{gaussian}.
}
  \item{alpha.level}{
	level, such that quantile \eqn{\tau=(1-}\code{alpha.level}\eqn{)/\gamma}. Default is 0.05.
}
  \item{M}{
	number of Monte Carlo Simulations to estimate the distribution \eqn{\Lambda}. Default is 1000.
}
  \item{qut.standardize}{
	standardize matrix X with a quantile-based standardization. Default is TRUE.
}
  \item{intercept}{
	should intercept(s) be fitted (default=TRUE) or set to zero (FALSE).
}
  \item{no.penalty}{
	unpenalized subset of covariates.
}
  \item{offset}{
	a vector of length \code{n} that is included in the linear predictor. Useful for the "poisson" family (e.g. log of exposure time), or for refining a model by starting at a current fit. Default is NULL.
}
  \item{bootstrap}{
set TRUE if it is desired to bootstrap matrix X when computing the Quantile Universal Threshold (Random scenario). Default is TRUE.
}
\item{beta0}{
	coefficients of the unpenalized covariates for generating the null data for the Quantile Universal Threshold. By default is \code{NA} and it is estimated using the unpenalized covatiriates and/or the intercept if TRUE.
	If it is desired to set \code{beta0} in advance, then it should be a vector of size the number of unpenalized covariates including the intercept if \code{intercept=TRUE}, in the same order. If there are not unpenalized covariates and \code{intercept=TRUE}, then it must be a real number.
}
\item{method}{
	Objective function for the zero thresholding. Select \code{lasso} for GLM-lasso or \code{sqrtlasso} for Square-root lasso.
}
}
\value{
  \item{lambda}{value of the Quantile Universal Threshold.}
  \item{Xnew}{standardized matrix X; \code{Xnew}=\code{X}\eqn{\times}\code{scale.factor}.}
  \item{scale.factor}{scale factor for \code{Xnew}.}
  \item{lambda.max}{smallest lambda that sets the lasso estimates to the zero vector.}
  \item{beta0 }{estimated value of the intercept when family is not \code{gaussian}.}
}
\references{
C. Giacobino, J. Diaz, S. Sardy, N. Hengartner. Quantile universal threshold for model selection. 2016
Jianqing Fan, Shaojun Guo and Ning Hao. Variance estimation using refitted cross-validation in ultrahigh dimensional regression. Journal of the Royal Statistical Society: Series B. 2012
Stephen Reid, Robert Tibshirani, and Jerome Friedman. A Study of Error Variance Estimation in Lasso Regression. 2013

}
\author{
Jairo Diaz
}


\seealso{
\code{\link{qut}}
}
\examples{
	X=matrix(rnorm(20*200),20,200)
	y=rnorm(20)+1
	lambda=lambdaqut(y,X,family=gaussian)
	
}