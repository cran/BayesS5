\name{result}
\alias{result}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Posterior inference results from the object of S5
}
\description{
Using the object of S5, the  maximum a posteriori (MAP) model, its posterior probability, and the marginal inclusion probabilities are provided.
}
\usage{
result(fit)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit}{
   an object of the 'S5' function.
}
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{hppm}{the MAP model }
   \item{hppm.prob}{the posterior probability of the MAP model}
    \item{marg.prob}{the marginal inclusion probabilities }
     \item{gam}{the binary vaiables of searched models by S5}
      \item{obj}{the corresponding log (unnormalized) posterior model probabilities}
      \item{post}{the corresponding (normalized) posterior model probabilities}
      \item{tuning}{the tuning parameter used in the model selection}
}
\references{
Shin, M., Bhattacharya, A., Johnson V. E. (2018) A Scalable Bayesian Variable Selection Using Nonlocal Prior Densities in Ultrahigh-dimensional Settings, Statistica Sinica,28(2): 1053-1078. 


Hans, C., Dobra, A., and West, M. (2007). Shotgun stochastic search for large p regression. Journal of the American Statistical Association, 102, 507-516.


Nikooienejad,A., Wang, W., and Johnson V.E. (2016). Bayesian variable selection for binary outcomes in high dimensional genomic studies using non-local priors. Bioinformatics, 32(9), 1338-45.

}
\author{
Shin Minsuk and Ruoxuan Tian
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\examples{
p=5000
n = 200

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
#str(res_pimom)
#print(res_pimom$hppm) 
#print(res_pimom$hppm.prob)  
#plot(res_pimom$marg.prob,ylim=c(0,1)) 
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
