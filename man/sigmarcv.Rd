\name{sigmarcv}
\alias{sigmarcv}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Variance estimation using refitted cross-validation
}
\description{
Variance estimation using refitted cross-validation in ultrahigh dimensional regression.
}
\usage{
sigmarcv(y, X, cv = FALSE, fit = NA, intercept = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
	response variable. Quantitative for family=\code{gaussian}, or family=\code{poisson} (non-negative counts). For family=\code{binomial} should be a factor with two levels.
}
  \item{X}{
	input matrix, of dimension n x p; each row is an observation vector.
}
  \item{cv}{
when FALSE, variance is estimated using Refitted Cross Validation in Fan et al. 2012; and when TRUE, it is estimated using cross validation as in Reid et al. 2013. Default is FALSE.
}
  \item{fit}{
A user supplied \code{glmnet} or {lars} object.  Typical usage is to leave it empty so that the program computes the regularization path using the algorithm selected in \code{type}. WARNING: use with care, if supplied, object options must match with user supplied options.
}
  \item{intercept}{
should intercept(s) be fitted (default=TRUE) or set to zero (FALSE).
}
}

\value{
	Estimator of \eqn{\sigma}
}
\references{
Variance estimation using refitted cross-validation in ultrahigh dimensional regression. Jianqing Fan, Shaojun Guo and Ning Hao. Journal of the Royal Statistical Society: Series B. 2012

}
\author{
Jianqing Fan, Shaojun Guo. Modified by Jairo Diaz.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
