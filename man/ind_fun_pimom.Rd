\name{ind_fun_pimom}
\alias{ind_fun_pimom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the log-marginal likelhood function based on piMoM priors
}
\description{
a log-marginal likelhood value of a model, based on the  piMoM prior on the regression coefficients and inverse gamma prior (0.01,0.01) on the variance.}
\usage{
ind_fun_pimom(X.ind,y,n,p,tuning)
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
Shin, M., Bhattacharya, A., Johnson V. E. (2018) A Scalable Bayesian Variable Selection Using Nonlocal Prior Densities in Ultrahigh-dimensional Settings, Statistica Sinica, 28(2), 1053-1078. 

  Johnson, V. E. and Rossell, D. (2012) Bayesian model selection in high-dimensional settings  , David, Journal of the American Statistical Association,  107 (498), 649-660.


}
\seealso{
\code{\link{ind_fun_g}}, \code{\link{ind_fun_pemom}}
}
