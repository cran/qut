\name{QUT-package}
\alias{QUT-package}
\docType{package}
\title{
Quantile Universal Threshold
}
\description{
Selection of a threshold parameter \eqn{\lambda} for GLM-lasso. The method consists in considering a null
model, finding the theoretical distribution of the threshold parameter under
the null, and setting \eqn{\lambda} to an upper quantile of that distribution. Although
this strategy does not use the data to select \eqn{\lambda} but simply considers
the behavior under the null model, it provides a theoretically and computationally
sound selection.
}
\details{
\tabular{ll}{
Package: \tab QUT\cr
Type: \tab Package\cr
License: \tab GPL-2\cr
}

}
\author{
Jairo Diaz <jairo.diaz@unige.ch>
}
\references{
Quantile universal threshold: model selection at
the detection edge for high-dimensional
generalized linear models, C. Giacobino, J. Diaz, S. Sardy, N. Hengartner.
}
\keyword{ QUT, GLM, threshold, lasso, regularization }


