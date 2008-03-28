\name{SimulateMixture}
\alias{SimulateMixture}

\title{Random Generation from a t Mixture Model} 

\description{
This function can be used to generate a sample from a multivariate \eqn{t} mixture model. 
}

\usage{
SimulateMixture(N, nu=100, mu, sigma, w)
}

\arguments{
\item{N}{The number of observations.}
\item{nu}{The degrees of freedom used for the \eqn{t} distribution.}
\item{mu}{A matrix of size \eqn{K \times P}{K x P}, where \eqn{K} is the number of clusters and \eqn{P} is the dimension, containing the \eqn{K} mean vectors.}
\item{sigma}{An array of dimension \eqn{K \times P \times P}{K x P x P}, containing the \eqn{K} covariance matrices.}
\item{w}{A vector of length \eqn{K}, containing the \eqn{K} cluster proportions.}
}

\value{
A matrix of size \eqn{N \times P}{N x P}.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\seealso{
\code{\link{flowClust}}
}

\examples{
% library(flowClust)
### Number of components
K <- 5
### Dimension
p <- 2
### Number of observations
n <- 200
Mu <- matrix(runif(K*p, 0, 20), K, p)
Sigma <- array(0, c(K, p, p))

for (k in 1:K)
{
    Sigma[k,,][outer(1:p, 1:p, ">")] <- runif(p*(p-1)/2,-.1,.1)
    diag(Sigma[k,,]) <- runif(p,0,1)
    ### Make sigma positive definite
    Sigma[k,,] <- Sigma[k,,] \%*\% t(Sigma[k,,])
}

### Generate the weights
w <- rgamma(K,10,1)
w <- w/sum(w)

y <- SimulateMixture(n, nu=4, Mu, Sigma, w)
}

\keyword{datagen}
