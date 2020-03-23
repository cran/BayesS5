\name{S5_additive}
\alias{S5_additive}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simplified shotgun stochastic search algorithm with screening (S5) for additive models
}
\description{
This is the Simplified Shotgun Stochastic Search with Screening (S5) for high-dimensonal Bayesian variable selection under nonparameteric additive models, which is considered in "Nonlocal Functional Priors for Nonparametric Hypothesis Testing and High-dimensional Model Selection" by Shin and Bhattacharya (2020+). This function utilizes the inverse moment nonlocal functional prior, and see Shin and Bhattacharya (2020+) for details.
}
\usage{
S5_additive(X, y, K=5, model, tuning, tem, ITER = 20, S = 20, C0 = 5, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
the covariate matrix (a standardization is recommneded for nonlocal priors).
}
  \item{y}{
a response variable. 
}
  \item{K}{
the degree of freedom for the B-spline basis
}

  \item{model}{
a model prior; Uniform or Bernoulli_Uniform. The default is Bernoulli_Uniform
}
  \item{tuning}{
a tuning parameter for the objective function (tau = g for the inverse moment prior).
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

Unlike "S5" function that requires an extra step to get more information of the computation procedure, this function provides full information of the results. 
}

\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
\item{GAM}{ the binary vaiables of searched models by S5}
\item{OBJ}{ the corresponding log (unnormalized) posterior probability}
\item{phi}{ the matrix of B-spline basis functions }
\item{Knots}{ the boundaries of knots used in generating the B-spline matrix }
\item{K}{ the degree of freedom of the B-spline basis.  }
\item{post}{the corresponding (normalized) posterior model probabilities}
\item{marg.prob}{the marginal inclusion probabilities }
\item{ind.MAP}{the selected variables from the MAP model}
\item{ind.MAP}{the selected variables whose marginal inclusion probability is larger than 0.5}
\item{hppm.prob}{the posterior probability of the MAP model}
\item{tuning}{the tuning parameter used in the model selection}
%% ...
}

\references{
Shin, M. and  Bhattacharya, A.(2020) Nonlocal Functional Priors for Nonparametric Hypothesis Testing and High-dimensional  Model Selection. 

Shin, M., Bhattacharya, A., Johnson V. E. (2018) A Scalable Bayesian Variable Selection Using Nonlocal Prior Densities in Ultrahigh-dimensional Settings, under revision in Statistica Sinica. 
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
p0 = 100
n0= 100

X = matrix(rnorm(n0*p0),n0,p0)
mu = 0.2*X[,1]^2 + sin(X[,2]*3) + cos(X[,3]*3) + 0.5*exp(X[,4])
y = mu + rnorm(n0)
X = scale(X)
y = as.vector(y)

### default setting
#fit_default = S5_additive(X,y)
#print(res_default$ind.hppm) # the MAP model 
#print(res_default$hppm.prob) # the posterior probability of the hppm 
#plot(res_default$marg.prob,ylim=c(0,1),ylab="marginal inclusion probability") 
# the marginal inclusion probability 
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{Simplified Shotgun Stochastic Search with Screening}% use one of  RShowDoc("KEYWORDS")
\keyword{Bayesian variable selection}% __ONLY ONE__ keyword per line
\keyword{nonlocal functional prior}
