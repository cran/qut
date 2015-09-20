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
lambdaqut(y, X, family = gaussian, alpha.level = "default", M = 1000, 
qut.standardize = TRUE, intercept = TRUE, no.penalty = NULL, offset = NULL,
bootstrap=TRUE)
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
	level, such that quantile \eqn{\tau=(1-}\code{alpha.level}\eqn{)/\gamma}. Default is \eqn{1/(\sqrt{\pi\log(p)})}.
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
}
\value{
  \item{lambda}{value of the Quantile Universal Threshold.}
  \item{Xnew}{standardized matrix X; \code{Xnew}=\code{X}\eqn{\times}\code{scale.factor}.}
  \item{scale.factor}{scale factor for \code{Xnew}.}
  \item{lambda.max}{smallest lambda that sets the lasso estimates to the zero vector.}
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