\name{predict.qut}
\alias{predict.qut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Make predictions from a "qut" object.
}
\description{
Similar to other predict methods, this function predicts fitted values from a fitted "qut" object
}
\usage{
\method{predict}{qut}(object, newx, mode = "glm", offset = NULL,...)
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
	x=matrix(rnorm(100*20),100,20)
	y1=x[,1]*10+rnorm(100)
	y2=as.numeric(binomial()$linkinv(x[,1]*10+rnorm(100))>0.5)

	fit1=qut(y1,x,family=gaussian,sigma=1)
	predict(fit1,newx=x[1:5,])
	predict(fit1,newx=x[1:5,],mode='lasso')
	y1[1:5]
	
	fit2=qut(y2,x,family=binomial)
	as.numeric(predict(fit2,newx=x[2:5,])>0.5)
	as.numeric(predict(fit2,newx=x[2:5,],mode='lasso')>0.5)
	y2[2:5]

}

