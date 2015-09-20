\name{sigmaqut}
\alias{sigmaqut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Estimation of \eqn{\sigma} based on the Quantile Universal Threshold
}
\description{
Estimation of \eqn{\sigma} using a two layer estimation scheme as in Refitted Cross Validation, by performing variable selection with the Quantile Universal Threshold, and obtaining the two estimations of sigma with the ordinary least squares estimator.
}
\usage{
sigmaqut(y, X, estimator = "unbiased", intercept = TRUE, 
alpha.level = "default", M = 1000, qut.standardize = TRUE, 
penalty.factor = rep(1, p), offset = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	response variable. Quantitative for family=\code{gaussian}, or family=\code{poisson} (non-negative counts). For family=\code{binomial} should be a factor with two levels.
}
  \item{X}{
     input matrix, of dimension n x p; each row is an observation vector.
}
  \item{estimator}{
	type of estimation of sigma when \code{sigma} = 'qut'. It can be equal to 'unbiased' (standard unbiased formula), or 'mle' (maximum likelihood formula).
}
  \item{intercept}{
	should intercept(s) be fitted (default=TRUE) or set to zero (FALSE).
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
  \item{penalty.factor}{
	separate penalty factors can be applied to each coefficient. As in \code{qut}.
}
  \item{offset}{
	a vector of length \code{n} that is included in the linear predictor. As in \code{qut}.
}
  \item{\dots}{
	other \code{glmnet} options.
}
}

\value{
Estimator of \eqn{\sigma}
}

\author{
Jairo Diaz
}
\note{
\code{\link{lambdaqut}},\code{\link{qut}}
}


