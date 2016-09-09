\name{slimsd}
\alias{slimsd}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Sparse Linear Regression using Nonsmooth Loss Functions and L1 Regularization
}
\description{
Function 'slim' of the 'flare package with the posibility of performing a different normalization
}
\usage{
slimsd(X, Y, lambda = NULL, nlambda = NULL, lambda.min.value = NULL, 
    lambda.min.ratio = NULL, rho = 1, method = "lq", q = 2, res.sd = FALSE, 
    prec = 1e-05, max.ite = 1e+05, verbose = FALSE,normalize=FALSE) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
	same as \code{slim}.
}
  \item{X}{
	same as \code{slim}.
}
  \item{lambda}{
same as \code{slim}.
}
  \item{nlambda}{
same as \code{slim}.
}
  \item{lambda.min.value}{
same as \code{slim}.
}
  \item{lambda.min.ratio}{
same as \code{slim}.
}
  \item{rho}{
same as \code{slim}.
}
  \item{method}{
same as \code{slim}.
}
  \item{q}{
same as \code{slim}.
}
  \item{res.sd}{
same as \code{slim}.
}
  \item{prec}{
same as \code{slim}.
}
  \item{max.ite}{
same as \code{slim}.
}
  \item{verbose}{
same as \code{slim}.
}
  \item{normalize}{
if TRUE, each variable is standardized as in \code{slim}, otherwise is left without standardization.
}
}

\value{
	same as \code{slim}.
}
\references{
Xingguo Li, Tuo Zhao, Xiaoming Yuan, Han Liu. The flare Package for High Dimensional Linear Regression and Precision Matrix Estimation in R. Journal of Machine Learning Research. 2015.
}
\author{
Xingguo Li, Tuo Zhao, Lie Wang, Xiaoming Yuan and Han Liu. Modified by Jairo Diaz.

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
