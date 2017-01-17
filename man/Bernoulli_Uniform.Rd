\name{Bernoulli_Uniform}
\alias{Bernoulli_Uniform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Bernoulli-Uniform model prior 
}
\description{
A mixture model prior with Bernoulli and uniform densities. See Scott and Berger (2010) for details. 
}
\usage{
Bernoulli_Uniform(ind,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ind}{
   an index set of variables in a model
}
\item{p}{
	the total number of covariates
}
}

\seealso{
	\code{\link{Uniform}}}  
%% ~~objects to See Also as \code{\link{help}}, ~~~


\references{
Scott, James G., and James O. Berger. "Bayes and empirical-Bayes multiplicity adjustment in the variable-selection problem." The Annals of Statistics 38.5 (2010): 2587-2619. }


\examples{
p = 5000
ind = 1:3 
m = Bernoulli_Uniform(ind,p)
print(m)
}




