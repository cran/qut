\name{qut}
\alias{qut}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit a low dimensional GLM with lasso using the Quantile Universal Threshold
}
\description{
Variable selection with GLM-lasso choosing the penalty parameter \eqn{\lambda} with the Quantile Universal Threshold. The procedure goes towards sparse estimation of the coefficients for good selection of the important predictors.
}
\usage{
qut(y,X,fit,family=gaussian,alpha.level='default',M=1000,qut.standardize=TRUE,
intercept=TRUE,offset=NULL,bootstrap=TRUE,sigma='qut',estimator='unbiased',
type=c('glmnet','lars'),lambda.seq=0,penalty.factor=rep(1,p),
lambda.min.ratio=ifelse(n<p,0.01,0.0001),nlambda=100,lambda=NULL,...)
}
\arguments{
  \item{y}{
	response variable. Quantitative for family=\code{gaussian}, or family=\code{poisson} (non-negative counts). For family=\code{binomial} should be a factor with two levels.
}
  \item{X}{
     input matrix, of dimension n x p; each row is an observation vector.
}
  \item{fit}{
     a user supplied \code{glmnet} or {lars} object.  Typical usage is to leave it empty so that the program computes the regularization path using the algorithm selected in \code{type}. WARNING: use with care, if supplied, object options must match with user supplied options.
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
  \item{offset}{
a vector of length \code{n} that is included in the linear predictor. Useful for the \code{poisson} family (e.g. log of exposure time), or for refining a model by starting at a current fit. Default is NULL.
}
  \item{bootstrap}{
set TRUE if it is desired to bootstrap matrix X when computing the Quantile Universal Threshold (Random scenario). Default is TRUE.
}
  \item{sigma}{
	standard deviation of the Gaussian errors. Used only if family=\code{gaussian}. When sigma = 'qut', it is estimated based on the Quantile Universal Threshold; when sigma = 'rcv', it is estimated using Refitted Cross Validation in Fan et al. 2012; and when sigma = 'cv', it is estimated using cross validation as in Reid et al. 2013. If sigma is a positive real number, then that value is used for the standard deviation. 
}
  \item{estimator}{
	type of estimation of sigma when \code{sigma} = 'qut'. It can be equal to 'unbiased' (standard unbiased formula), or 'mle' (maximum likelihood formula).
}
  \item{type}{
algorithm for constructing the regularization path. It can be \code{lars} (type='lars') or \code{glmnet} (type='glmnet'). If family is not gaussian, penalty.factor is different from default, or offset different from NULL, glmnet will be always used. Default is 'glmnet'.
}
  \item{lambda.seq}{
preset lambda sequence when type = 'glmnet'. If \code{lambda.seq}<2 the sequence of lambdas goes decreasing from \code{lambda.max} to \code{lambda.qut}. if \code{lambda.seq}= 0 or 1, \code{lambda} sequence is equispaced on the log-scale, with \code{nlambda}=100 for 0 and \code{nlambda}=n for 1. Use \code{lambda.seq}=2 for \code{glmnet} default options. Default is 0.
}
  \item{penalty.factor}{
separate penalty factors can be applied to each coefficient. This is a number that multiplies lambda to allow differential shrinkage. Can be 0 for some variables, which implies no shrinkage, and that variable is always included in the model. Default is 1 for all variables (and implicitly infinity for variables listed in exclude). Note: the penalty factors are internally rescaled to sum to n, and the lambda sequence will reflect this change.
}
  \item{lambda.min.ratio}{
smallest value for lambda, as a fraction of \code{lambda.max}. As in \code{glmnet}.
}
  \item{nlambda}{
the number of \code{lambda}. As in \code{glmnet}. Default is 100.
}
  \item{lambda}{
	a user supplied \code{lambda} sequence. As in \code{glmnet}.
}

  \item{\dots}{
	\code{glmnet} or \code{lars} options.
}
}

\value{
  \item{lambda}{value of the Quantile Universal Threshold.}
  \item{fit }{object fitted by \code{glmnet} or \code{lars}.}
  \item{beta }{coefficients obtained with the Quantile Universal Threshold.}
  \item{betaglm }{coefficients obtained fitting GLM with the non zero coefficients in \code{beta}.}
  \item{family}{response type}
  \item{sigma}{standard deviation estimate of the errors (when family=\code{gaussian})}
  \item{scale.factor}{scale factor used for standardizing \eqn{X}.}
}
\references{
Quantile universal threshold selection with an
application in generalized linear models with LASSO, C. Giacobino, J. Diaz, S. Sardy, N. Hengartner.

Quantile universal threshold: model selection at the detection edge for high-dimensional linear regression, J. Diaz, S. Sardy.

Variance estimation using refitted cross-validation in ultrahigh dimensional regression. Jianqing Fan, Shaojun Guo and Ning Hao. Journal of the Royal Statistical Society: Series B. 2012

A Study of Error Variance Estimation in Lasso Regression. Stephen Reid, Robert Tibshirani, and Jerome Friedman. 2013

}
\author{
Jairo Diaz
}


\seealso{
 \code{\link{lambdaqut}}
}
\examples{

X=matrix(rnorm(50*500),50,500)
beta=c(rep(10,5),rep(0,500-5))
y=X \%*\% beta+rnorm(50)

outqut=qut(y,X,type='glmnet',family=gaussian,sigma=1) #Fitting with qut
betaqut=outqut$beta[-1]

outcv=cv.glmnet(X,y,family='gaussian') #fitting with Cross-Validation
betacv=coef(outcv$glmnet.fit,s=outcv$lambda.min)[-1]

results=rbind( 	c(sum(betaqut[1:5]!=0),sum(betaqut[-(1:5)]!=0)),
				c(sum( betacv[1:5]!=0), sum(betacv[-(1:5)]!=0)) )
colnames(results)=c('True Positive','False Positive')
rownames(results)=c('qut','cv')

print(results)

}