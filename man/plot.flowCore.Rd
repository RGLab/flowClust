\name{plot,flowFrame,tmixFilterResult-method}
\docType{methods}
\alias{plot,flowFrame,tmixFilterResult-method}
\alias{plot,flowFrame,tmixFilterResultList-method}
\alias{plot,tmixFilterResult-method}
\alias{plot.flowFrame.tmixFilterResult}
\alias{plot.flowFrame}
\alias{plot.tmixFilterResult}

\title{Scatterplot / 1-D Density Plot of Filtering (Clustering) Results} 

\description{
Depending on the dimensions specified, this method generates either a scatterplot or a one-dimensional density plot (histogram) based on the robust model-based clustering results.
}

\usage{
\S4method{plot}{flowFrame,tmixFilterResult}(x, y, z=NULL, \dots)
}

\arguments{
\item{x}{Object of class \code{flowFrame}.  This is the data object on which \code{\link[=tmixFilter]{filter}} was performed.}
\item{y}{Object of class \code{tmixFilterResult} or \code{tmixFilterResultList} returned from running \code{\link[=tmixFilter]{filter}}.}
\item{z}{A character vector of length one or two containing the name(s) of the variable(s) selected for the plot.  If it is of length two, a scatterplot will be generated.  If it is of length one, a 1-D density plot will be made.  If it is unspecified, the first one/two variable(s) listed in \code{y@varNames} will be used.}
\item{\dots}{All optional arguments passed to the \code{\link[=plot.flowClust]{plot}} or \code{\link[=hist.flowClust]{hist}} method with signature \code{'flowClust'}.  Note that arguments \code{x}, \code{data} and \code{subset} have already been provided by \code{y}, \code{x} and \code{z} above respectively.}
}

\note{
This \code{plot} method is designed such that it resembles the argument list of the \code{plot} method defined in the \pkg{flowCore} package.  The actual implementation is done through the \code{\link[=plot.flowClust]{plot}} or \code{\link[=hist.flowClust]{hist}} method defined for a \code{flowClust} object.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link[=tmixFilter]{filter}}, \code{\link[=plot.flowClust]{plot}}, \code{\link[=hist.flowClust]{hist}}
}

\keyword{graphs}





