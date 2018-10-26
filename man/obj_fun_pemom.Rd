\name{obj_fun_pemom}
\alias{obj_fun_pemom}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the log posterior distribution based on peMoM priors and inverse gamma prior (0.01,0.01)
}
\description{
a log posterior density value at regression coefficients of a model, based on the  peMoM prior on the regression coefficients  and inverse gamma prior (0.01,0.01) on the variance.}
\usage{
obj_fun_pemom(ind,X,y,n,p,tuning)
}
\arguments{
    \item{ind}{
      the index set of  a model
      }
      
   \item{X}{
   	the covariates
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


 Rossell, D., Telesca, D., and Johnson, V. E. (2013) High-dimensional Bayesian classifiers using non-local priors, Statistical Models for Data Analysis, 305-313.
}

\seealso{
\code{\link{obj_fun_g}}, \code{\link{obj_fun_pimom}}
}
