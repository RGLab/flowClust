\name{posterior}
\alias{posterior}
\alias{importance}
\alias{uncertainty}
\alias{getEstimates}

\title{Various Functions for Retrieving Information from Clustering Results} 

\description{
Various functions are available to retrieve the posterior probabilities of clustering memberships \eqn{z} (\code{posterior}), the \dQuote{weights} \eqn{u} (\code{importance}), the uncertainty (\code{uncertainty}), and the estimates of the cluster proportions, means and variances (\code{getEstimates}) resulted from the clustering (filtering) operation.
}

\usage{
posterior(object, assign=FALSE)
importance(object, assign=FALSE)
uncertainty(object)
getEstimates(object, data)
}

\arguments{
\item{object}{Object returned from \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}}.}
\item{assign}{A logical value.  If \code{TRUE}, only the quantity (\code{z} for \code{posterior} or \code{u} for \code{importance}) associated with the cluster to which an observation is assigned will be returned.  Default is \code{FALSE}, meaning that the quantities associated with all the clusters will be returned.}
\item{data}{A numeric vector, matrix, data frame of observations, or object of class \code{flowFrame}; an optional argument.  This is the object on which \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}} was performed.}
}

\details{
These functions are written to retrieve various slots contained in the object returned from the clustering operation.  \code{posterior} and \code{importance} provide a means to conveniently retrieve information stored in \code{object@z} and \code{object@u} respectively.  \code{uncertainty} is to retrieve \code{object@uncertainty}. \code{getEstimates} is to retrieve information stored in \code{object@mu} (transformed back to the original scale) and \code{object@w}; when the data object is provided, an approximate variance estimate (on the original scale, obtained by performing one M-step of the EM algorithm without taking the Box-Cox transformation) will also be computed.
}

\value{
Denote by \eqn{K} the number of clusters, \eqn{N} the number of observations, and \eqn{P} the number of variables.  For \code{posterior} and \code{importance}, a matrix of size \eqn{N \times K}{N x K} is returned if \code{assign=FALSE} (default).  Otherwise, a vector of size \eqn{N} is outputted.  \code{uncertainty} always outputs a vector of size \eqn{N}.  \code{getEstimates} returns a list with named elements, \code{proportions}, \code{locations} and, if the data object is provided, \code{dispersion}.  \code{proportions} is a vector of size \eqn{P} and contains the estimates of the \eqn{K} cluster proportions.  \code{locations} is a matrix of size \eqn{K \times P}{K x P} and contains the estimates of the \eqn{K} mean vectors transformed back to the original scale (i.e., \code{rbox(object@mu, object@lambda)}).  \code{dispersion} is an array of dimensions \eqn{K \times P \times P}{K x P x P}, containing the approximate estimates of the \eqn{K} covariance matrices on the original scale.
}

% \note{
% If the \code{assign} argument is set to \code{TRUE}, only the quantities corresponding to assigned observations will be returned.  Quantities corresponding to unassigned observations (outliers and filtered observations) will be reported as \code{NA}.  Hence, A change in the rule to call outliers will incur a change in the number of \code{NA} values returned.
% }

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}}, \code{\link{Map}}
}

\examples{
res <- flowClust(iris[,1:4], K=3)
posterior(res)
posterior(res, assign=TRUE)
importance(res)
importance(res, assign=TRUE)
uncertainty(res)
getEstimates(res)
}

\keyword{cluster}
