\name{ind_fun_pemom}
\alias{ind_fun_pemom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the log-marginal likelhood function based on peMoM priors and inverse gamma prior (0.01,0.01)
}
\description{
a log-marginal likelhood value of a model, based on the  peMoM prior on the regression coefficients  and inverse gamma prior (0.01,0.01) on the variance.}
\usage{
ind_fun_pemom(X.ind,y,n,p,tuning)
}
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
Shin, M., Bhattacharya, A., Johnson V. E. (2017+) A Scalable Bayesian Variable Selection Using Nonlocal Prior Densities in Ultrahigh-dimensional Settings, accepted in Statistica Sinica. 

 Rossell, D., Telesca, D., and Johnson, V. E. (2013) High-dimensional Bayesian classifiers using non-local priors, Statistical Models for Data Analysis, 305-313.
}

\seealso{
\code{\link{ind_fun_g}}, \code{\link{ind_fun_pimom}}
}
