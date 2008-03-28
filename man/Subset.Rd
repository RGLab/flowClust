\name{Subset,flowClust-method}
\docType{methods}
\alias{Subset,ANY,flowClust-method}
\alias{Subset,flowClust-method}
\alias{Subset.flowClust}
\alias{Subset,flowFrame,tmixFilterResult-method}
\alias{Subset.flowFrame}
\alias{Subset.tmixFilterResult}
\alias{Subset.flowFrame.tmixFilterResult}
\alias{Subset}

\title{Subsetting Data Based on Clustering Results} 

\description{
This method returns a subset of data upon the removal of outliers identified from the clustering (filtering) operations.
}

\usage{
\S4method{Subset}{ANY,flowClust}(x, subset, \dots)

\S4method{Subset}{flowFrame,tmixFilterResult}(x, subset, \dots)
}

\arguments{
\item{x}{A numeric vector, matrix, data frame of observations, or object of class \code{flowFrame}.  This is the object on which \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}} was performed.}
\item{subset}{Object returned from \code{flowClust} or \code{filter}.}
\item{\dots}{Further arguments to be passed to or from other methods.  The selection of columns can be done by specifying the optional \code{select} argument that takes a numeric (not supported when \code{x} is of class \code{flowFrame}) or character vector.}
}

\value{
An object which is a subset of \code{x}.  It also retains the same class as \code{x}.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{split}}, \code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}}
}

\keyword{manip}
