\name{result_est_MAP}
\alias{result_est_MAP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior inference results from the object of S5
}
\description{
Using the object of S5, the maximum a posteriori (MAP) estimator  and Bayesian Model Averaged (BMA) estimators of the regression coefficients are provided.
}
\usage{
result_est_MAP(res,X,y,obj_fun,verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{res}{
   an object of the 'S5' function.
}
\item{X}{the covariates.}
\item{y}{the response varaible.}
\item{obj_fun}{the negative log (unnormalized) posterior density when a model is given.}
\item{verbose}{logical; default is TRUE.}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{intercept.MAP}{the MAP estimator of the intercept. }
   \item{beta.MAP}{the MAP estimator of the regression coefficients.}
  \item{sig.MAP}{the MAP estimator of the regression variance.}
  \item{intercept.BMA}{the Baeysian model averaged estimator of the intercept. }
   \item{beta.BMA}{the Bayesian model averaged estimator of the regression coefficients.}
  
  }
\references{
Shin, M., Bhattacharya, A., Johnson V. E. (2018) A Scalable Bayesian Variable Selection Using Nonlocal Prior Densities in Ultrahigh-dimensional Settings, Statistica Sinica. 

Hans, C., Dobra, A., and West, M. (2007). Shotgun stochastic search for large p regression. Journal of the American Statistical Association, 102, 507-516.


Nikooienejad,A., Wang, W., and Johnson V.E. (2016). Bayesian variable selection for binary outcomes in high dimensional genomic studies using non-local priors. Bioinformatics, 32(9), 1338-45.

}
\author{
Shin Minsuk and Ruoxuan Tian
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
p=5000
n = 100

indx.beta = 1:5
xd0 = rep(0,p);xd0[indx.beta]=1
bt0 = rep(0,p); 
bt0[1:5]=c(1,1.25,1.5,1.75,2)*sample(c(1,-1),5,replace=TRUE)
xd=xd0
bt=bt0
X = matrix(rnorm(n*p),n,p)
y = X\%*\%bt0 + rnorm(n)*sqrt(1.5)
X = scale(X)
y = y-mean(y)
y = as.vector(y)

### piMoM  
#C0 = 2 # the number of repetitions of S5 algorithms to explore the model space
#tuning = 10 # tuning parameter
#tuning = hyper_par(type="pimom",X,y,thre = p^-0.5)
#print(tuning)
#ind_fun = ind_fun_pimom # choose the prior on the regression coefficients (pimom in this case)
#model = Bernoulli_Uniform # choose the model prior 
#tem =  seq(0.4,1,length.out=20)^2 # the sequence of the temperatures

#fit_pimom = S5(X,y,ind_fun=ind_fun,model = model,tuning=tuning,tem=tem,C0=C0)
#fit_pimom$GAM # the searched models by S5
#fit_pimom$OBJ # the corresponding log (unnormalized) posterior probability

#res_pimom = result(fit_pimom)
#est.MAP = result_est_MAP(res_pimom,X,y,obj_fun_pimom,verbose=TRUE)
#plot(est.MAP$beta.MAP,est.MAP$beta.BMA)
#abline(0,1,col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Simplified Shotgun Stochastic Search with Screening}% use one of  RShowDoc("KEYWORDS")
\keyword{Shotgun Stochastic Search}% use one of  RShowDoc("KEYWORDS")
\keyword{Bayesian variable selection}% __ONLY ONE__ keyword per line
\keyword{g-prior}
\keyword{nonlocal prior}
\keyword{piMoM}
\keyword{peMoM}
