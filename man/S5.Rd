\name{S5}
\alias{S5}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simplified shotgun stochastic search algorithm with screening (S5)
}
\description{
The Simplified Shotgun Stochastic Search with Screening (S5) is proposed by Shin et al (2018), which is a scalable stochastic search algorithm for high-dimensonal Bayesian variable selection.  It is a modified version of the Shotgun Stochasitic Search (SSS, Hans et al., 2007), aimed at rapidly identifying regions of high posterior probability and finding the maximum a posteriori (MAP) model. Also, the S5 provides an approximation of posterior probability of each model (including the marginal inculsion probabilities). For details, see Shin et al. (2018)
}
\usage{
S5(X, y, ind_fun, model, tuning, tem, ITER = 20, S = 20, C0 = 5, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the covariate matrix (a standardization is recommneded for nonlocal priors).
}
  \item{y}{
a response variable. 
}
  \item{ind_fun}{
  a log-marginal likelihood function of models, which is resulted from a pred-specified priors on the regression coefficients. The default is "piMoM". See the example below for details.
}
  \item{model}{
a model prior; Uniform or Bernoulli_Uniform. The default is Bernoulli_Uniform
}
  \item{tuning}{
a tuning parameter for the objective function (tau for piMoM and peMoM priors; g for the g-prior).
}
  \item{tem}{
a temperature schedule. The default is seq(0.4,1,length.out=20)^-2.
}
  \item{ITER}{
 the number of iterations in each temperature; default is 20.
}
  \item{S}{
a screening size of variables; default is 20.
}
  \item{C0}{
a number of repetition of  the S5 algorithm C0 times,default is 2. When the total number of variables is huge and real data sets are considered, using a large number of C0 is recommended, e.g., C0=5.  
}
 \item{verbose}{
 	if TRUE, the function prints the currnet status of the S5 in each temperature; the default is TRUE.
}
 }
 
\details{
Using the S5 (Shin et al., 2018), you will get all the models searched by S5 algorithm, and their corresponding log (unnormalized) posterior probabilities, and also this function can receive searched model for g-prior,piMoM,and peMoM. 

After obtaining the object of the S5 function, by using the 'result' function, you can obtain the posterior probabilities of the searched models including the MAP model and  the marginal inclusion probabilities of each variable.

By using the procedure of Nikooienejad et al. (2016), the 'hyper_par' function chooses the tuning parameter for nonlocal priors (piMoM or peMoM priors).
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{GAM}{ the binary vaiables of searched models by S5}
\item{OBJ}{ the corresponding log (unnormalized) posterior probability}
\item{tuning}{the tuning parameter used in the model selection}
%% ...
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

\seealso{
	\code{\link{result}}, \code{\link{S5_parallel}},  \code{\link{SSS}}  
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
p0 = 5000
n0= 100

indx.beta = 1:5
xd0 = rep(0,p0);xd0[indx.beta]=1
bt0 = rep(0,p0); 
bt0[1:5]=c(1,1.25,1.5,1.75,2)*sample(c(1,-1),5,replace=TRUE)
xd=xd0
bt=bt0
X = matrix(rnorm(n0*p0),n0,p0)
y = crossprod(t(X),bt0) + rnorm(n0)*sqrt(1.5)
X = scale(X)
y = y-mean(y)
y = as.vector(y)

### default setting
#fit_default = S5(X,y)
#res_default = result(fit_default)
#print(res_default$hppm) # the MAP model 
#print(res_default$hppm.prob) # the posterior probability of the hppm 
#plot(res_default$marg.prob,ylim=c(0,1),ylab="marginal inclusion probability") 
# the marginal inclusion probability 

### Nonlocal prior (piMoM prior) by S5
#C0 = 1 # the number of repetitions of S5 algorithms to explore the model space
#tuning = hyper_par(type="pimom",X,y,thre = p^-0.5)  
# tuning parameter selection for nonlocal priors
#print(tuning) 

#ind_fun = ind_fun_pimom # the log-marginal likelihood of models based on piMoM prior
#model = Bernoulli_Uniform 
# the log-marginal likelihood of models based on piMoM prior 
#tem =  seq(0.4,1,length.out=20)^2 
# the temperatures schedule
#fit_pimom = S5(X,y,ind_fun=ind_fun,model=model,tuning=tuning,tem=tem,C0=C0)


#fit_pimom$GAM # the searched models by S5
#fit_pimom$OBJ # the corresponding log (unnormalized) posterior probability

#res_pimom = result(fit_pimom)
#str(res_pimom)
#print(res_pimom$hppm) # the MAP model 
#print(res_pimom$hppm.prob) 
# the posterior probability of the hppm 
#plot(res_pimom$marg.prob,ylim=c(0,1),ylab="marginal inclusion probability") 
# the marginal inclusion probability 


### Get the estimated regression coefficients from Bayesian Model Avaeraging (BMA)
#est.LS = result_est_LS(res_pimom,X,y) # Averged over the Least Square estimators of the models.
#est.MAP = result_est_MAP(res_pimom,X,y,obj_fun_pimom,verbose=TRUE) 
# Averged over the maximum posteriori (MAP) estimators of the models.
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{Simplified Shotgun Stochastic Search with Screening}% use one of  RShowDoc("KEYWORDS")
\keyword{Shotgun Stochastic Search}% use one of  RShowDoc("KEYWORDS")
\keyword{Bayesian variable selection}% __ONLY ONE__ keyword per line
\keyword{g-prior}
\keyword{nonlocal prior}
\keyword{piMoM}
\keyword{peMoM}
