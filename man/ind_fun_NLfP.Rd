\name{ind_fun_NLfP}
\alias{ind_fun_NLfP}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
the log-marginal likelhood function based on the invers moment functional priors and inverse gamma prior (0.01,0.01)
}
\description{
a log-marginal likelhood value of a model, based on the  peMoM prior on the regression coefficients  and inverse gamma prior (0.01,0.01) on the variance.}
\usage{
ind_fun_NLfP(ind2, y, phi, n, p, K, IP.phi, C.prior1, tuning)
}
\arguments{
    \item{ind2}{
      the index of covariates in a model
      }
   \item{y}{
   	the response variable
   }
   \item{phi}{
   	the B-spline basis
   }
   \item{n}{
   	the sample size   
   	}   
   \item{p}{
   	the total number of covariates
   }
   \item{K}{
   	the degree of freedom for the B-spline basis
   }
   \item{IP.phi}{
   	the projection matrix on the null space; Q matrix in Shin and Bhattacharya (2020+)
   }
   \item{C.prior1}{
   	the logarithm of the normalizing constant of the nonlocal funcitonal prior
   }
   \item{tuning}{
   	a value of the tuning parameter
   }

}

\references{
Shin, M. and  Bhattacharya, A.(2020) Nonlocal Functional Priors for Nonparametric Hypothesis Testing and High-dimensional  Model Selection. 
}
