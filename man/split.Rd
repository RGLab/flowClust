\name{split,flowClust-method}
\docType{methods}
\alias{split,flowClust-method}
\alias{split.flowClust}
\alias{split,data.frame,flowClust-method}
\alias{split,matrix,flowClust-method}
\alias{split,vector,flowClust-method}
\alias{split,flowFrame,flowClust-method}
\alias{split,flowFrame,tmixFilterResult-method}
\alias{split,data.frame,flowClustList-method}
\alias{split,matrix,flowClustList-method}
\alias{split,vector,flowClustList-method}
\alias{split,flowFrame,flowClustList-method}
\alias{split,flowFrame,tmixFilterResultList-method}
\alias{split.flowFrame}
\alias{split.tmixFilterResult}
\alias{split}

\title{Splitting Data Based on Clustering Results} 

\description{
This method splits data according to results of the clustering (filtering) operation.  Outliers identified will be removed by default.
}

\section{Usage}{
split(x, f, drop=FALSE, population=NULL, split=NULL, rm.outliers=TRUE, \dots)
}

\section{Arguments}{
\item{x}{A numeric vector, matrix, data frame of observations, or object of class \code{flowFrame}.  This is the object on which \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}} was performed.}
\item{f}{Object returned from \code{flowClust} or \code{filter}.}
\item{drop}{A logical value indicating whether to coerce a column matrix into a vector, if applicable.  Default is \code{FALSE}, meaning that a single-column matrix will be retained.}
\item{population}{An optional argument which specifies how to split the data.  If specified, it takes a list object with named or unnamed elements each of which is a numeric vector specifying which clusters are included.  If this argument is left unspecified, the data object will be split into \code{K} subsets each of which is formed by one out of the \code{K} clusters used to model the data.  See examples for more details.}
\item{split}{This argument is deprecated.  Should use \code{population} instead.}
\item{rm.outliers}{A logical value indicating whether outliers are removed or not.}
\item{\dots}{Further arguments to be passed to or from other methods.}
}

\value{
A list object with elements each of which is a subset of \code{x} and also retains the same class as \code{x}.  If the \code{split} argument is specified with a list of named elements, those names will be used to name the corresponding elements in the resultant list object.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{Subset}}, \code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}}
}

\keyword{manip}
