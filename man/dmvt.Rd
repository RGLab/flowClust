\name{dmvt}
\alias{dmvt}

\title{Density of the Multivariate t Distribution with Box-Cox Tranformation} 

\description{
This function computes the densities at the inputted points of the multivariate \eqn{t} distribution with Box-Cox transformation.
}

\usage{
dmvt(x, mu, sigma, nu, lambda, log=FALSE)
}

\arguments{
\item{x}{A matrix or data frame of size \eqn{N \times P}{N x P}, where \eqn{N} is the number of observations and \eqn{P} is the dimension.  Each row corresponds to one observation.}
\item{mu}{A numeric vector of length \eqn{P} specifying the mean.}
\item{sigma}{A matrix of size \eqn{P \times P}{P x P} specifying the covariance matrix.}
\item{nu}{The degrees of freedom used for the \eqn{t} distribution.  If \code{nu=Inf}, Gaussian distribution will be used.}
\item{lambda}{The Box-Cox transformation parameter.  If missing, the conventional \eqn{t} distribution without transformation will be used.}
\item{log}{A logical value.  If \code{TRUE} then the logarithm of the densities is returned.}
}

\value{
A list with the following components:
\item{value}{A vector of length \eqn{N} containing the density values.}
\item{md}{A vector of length \eqn{N} containing the Mahalanobis distances.}
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

% \examples{
% 
% }

\keyword{distribution}





