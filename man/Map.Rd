\name{Map,flowClust-method}
\docType{methods}
\alias{Map,flowClust-method}
\alias{Map,flowClustList-method}
\alias{Map.flowClust}
\alias{Map}

\title{Cluster Assignment Based on Clustering Results} 

\description{
This method performs cluster assignment according to the posterior probabilities of clustering memberships resulted from the clustering (filtering) operations.  Outliers identified will be left unassigned by default.
}

\usage{
\S4method{Map}{flowClust}(f, rm.outliers=TRUE, \dots)
}

\arguments{
\item{f}{Object returned from \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}}.}
\item{rm.outliers}{A logical value indicating whether outliers will be left unassigned or not.}
\item{\dots}{Further arguments to be passed to or from other methods.}
}

\value{
A numeric vector of size \eqn{N} (the number of observations) indicating to which cluster each observation is assigned.  Unassigned observations will be labelled as \code{NA}.
}

\note{
Even if \code{rm.outliers} is set to \code{FALSE}, \code{NA} may still appear in the resultant vector due to the filtered observations; see the descriptions about the \code{min.count}, \code{max.count}, \code{min} and \code{max} arguments of \code{\link{flowClust}}.
}


\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}}, \code{\link{posterior}}
}

\examples{
data(rituximab)

### cluster the data using FSC.H and SSC.H
res1 <- flowClust(rituximab, varNames=c("FSC.H", "SSC.H"), K=1)

### remove outliers before proceeding to the second stage
rituximab2 <- Subset(rituximab, res1)

### cluster the data using FL1.H and FL3.H (with 3 clusters)
res2 <- flowClust(rituximab2, varNames=c("FL1.H", "FL3.H"), K=3)
summary(res2)

# to show the cluster assignment of observations
table(Map(res2))
}

\keyword{cluster}
