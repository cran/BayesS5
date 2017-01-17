\name{ind_fun_g}
\alias{ind_fun_g}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Zellner's g-prior 
}
\description{
a log-marginal likelhood value of a model, based on the  Zellner's g-prior on the regression coefficients.
}
\usage{
ind_fun_g(X.ind,y,n,p,tuning)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X.ind}{
      the subset of covariates in a model
      }
   \item{y}{
   	the response variable
   }   
   \item{n}{
   	the sample size   
   	}   
   \item{p}{
   	the total number of covariates
   }   
   \item{tuning}{
   	a value of the tuning parameter
   }

}


\references{
	Zellner, Arnold. "On assessing prior distributions and Bayesian regression analysis with g-prior distributions." Bayesian inference and decision techniques: Essays in Honor of Bruno De Finetti 6 (1986): 233-243.
}

\author{
Shin Minsuk and Ruoxuan Tian
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{ind_fun_pimom}}, \code{\link{ind_fun_g}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}

\examples{
#p=5000
p = 10
n = 200

indx.beta = 1:5
xd0 = rep(0,p);xd0[indx.beta]=1
bt0 = rep(0,p); 
bt0[1:5]=c(1,1.25,1.5,1.75,2)*sample(c(1,-1),5,replace=TRUE)
xd=xd0
bt=bt0
X = matrix(rnorm(n*p),n,p)
y = crossprod(t(X),bt0) + rnorm(n)*sqrt(1.5)
X = scale(X)
y = y-mean(y)
y = as.vector(y)

C0 = 1 # the number of repetitions of S5 algorithms to explore the model space
tuning = p^2 # tuning parameter g for g-prior
ind_fun = ind_fun_g # choose the pror on the regression coefficients (g-prior in this case)
model = Uniform #choose the model prior (Uniform prior in this cases)
tem =  seq(0.4,1,length.out=20)^2 # the sequence of the temperatures

fit_g = S5(X,y,ind_fun=ind_fun,model=model, tuning=tuning,tem=tem,C0=C0)
}
\keyword{Shotgun Stochastic Search}% use one of  RShowDoc("KEYWORDS")
\keyword{Bayesian variable selection}% __ONLY ONE__ keyword per line
\keyword{g-prior}
\keyword{nonlocal prior}
\keyword{piMoM}
\keyword{peMoM}
