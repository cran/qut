\name{predict.qut}
\alias{predict.qut}
\alias{coef.qut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make predictions from a "qut" object.
}
\description{
Similar to other predict methods, this function predicts fitted values from a fitted "qut" object
}
\usage{
\method{predict}{qut}(object, newx, mode = "glm", offset = NULL,...)
\method{coef}{qut}(object, mode = "glm",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
	fitted "qut" model object.
}
  \item{newx}{
	matrix of new values for X at which predictions are to be made. Must be a matrix.
}
  \item{mode}{
	make predictions with lasso coefficients (type=\code{lasso}) or with fitted glm coefficients (type=\code{glm}). Default is \code{glm}.
}
  \item{offset}{
	if an offset is used in the fit, then one must be supplied for making predictions
}
  \item{...}{
	not used. Other arguments to predict.
}
}

\value{
	a vector/matrix of fitted values
}

\author{
	Jairo Diaz
}

\seealso{
	\code{\link{qut}}
}
\examples{
	set.seed(1234)
	x=matrix(rnorm(200*20),200,20)
	y1=x[,1]*10+rnorm(100)

	fit1=qut(y1,x,family=gaussian,sigma=1)
	predict(fit1,newx=x[1:5,])
	predict(fit1,newx=x[1:5,],mode='lasso')
	y1[1:5]
	
	coef(fit1,mode='lasso')
	coef(fit1,mode='glm')

}

