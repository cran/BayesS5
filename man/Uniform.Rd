\name{Uniform}
\alias{Uniform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Uniform model prior
}
\description{
A uniform model prior that assigns the same prior mass on each model.
}
\usage{
Uniform(ind,p)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ind}{
    the index set of variables in a model
}

\item{p}{
	the total number of covariates
}
}

\examples{
ind = 1:3 
m = Uniform(ind,p)
print(m)
}


