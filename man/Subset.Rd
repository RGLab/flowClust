\name{Subset,flowClust-method}
\docType{methods}
\alias{Subset,flowClust-method}
\alias{Subset.flowClust}
\alias{Subset.flowFrame}
\alias{Subset.tmixFilterResult}
\alias{Subset.flowFrame.tmixFilterResult}
\alias{Subset,flowFrame,flowClust-method}
\alias{Subset,flowFrame,tmixFilterResult-method}
\alias{Subset,data.frame,flowClust-method}
\alias{Subset,matrix,flowClust-method}
\alias{Subset,vector,flowClust-method}
\alias{Subset,ANY,flowClustList-method}
\alias{Subset,flowFrame,tmixFilterResultList-method}
\alias{Subset}

\title{Subsetting Data Based on Clustering Results} 

\description{
This method returns a subset of data upon the removal of outliers identified from the clustering (filtering) operations.
}

\section{Usage}{
Subset(x, subset, \dots)
}

\arguments{
\item{x}{A numeric vector, matrix, data frame of observations, or object of class \code{flowFrame}.  This is the object on which \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}} was performed.}
\item{subset}{Object returned from \code{flowClust} or \code{filter}.}
\item{\dots}{Further arguments to be passed to or from other methods.}
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

