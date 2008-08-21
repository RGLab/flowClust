\name{hist,flowClust-method}
\docType{methods}
\alias{hist,flowClust-method}
\alias{hist.flowClust}

\title{1-D Density Plot (Histogram) of Clustering Results}

\description{
This method generates a one-dimensional density plot for the specified dimension (variable) based on the robust model-based clustering results.  A histogram of the actual data or cluster assignment is optional for display.
}

\usage{
\S4method{hist}{flowClust}(x, data=NULL, subset=1, include=1:(x@K), histogram=TRUE, 
     labels=TRUE, xlim=NULL, ylim=NULL, xlab=NULL, ylab="Density", 
     main=NULL, breaks=50, col=NULL, pch=20, cex=0.6, \dots)
}

\arguments{
\item{x}{Object returned from \code{\link{flowClust}} or from running \code{\link[=filter.flowFrame]{filter}} on a \code{flowFrame} object.}
\item{data}{A numeric vector, matrix, data frame of observations, or object of class \code{flowFrame}. This is the object on which \code{flowClust} or \code{filter} was performed.}
\item{subset}{An integer indicating which variable is selected for the plot.  Alternatively, a character string containing the name of the variable is allowed if \code{x@varNames} is not \code{NULL}.}
\item{include}{A numeric vector specifying which clusters are shown on the plot.  By default, all clusters are included.}
\item{histogram}{A logical value indicating whether a histogram of the actual data is made in addition to the density plot or not.}
\item{labels}{A logical value indicating whether information about cluster assignment is shown or not.}
\item{xlim}{The range of \eqn{x}-values for the plot.  If \code{NULL}, the data range will be used.}
\item{ylim}{The range of \eqn{y}-values for the plot.  If \code{NULL}, an optimal range will be determined automatically.}
\item{xlab, ylab}{Labels for the \eqn{x}- and \eqn{y}-axes respectively.}
\item{main}{Title of the plot.}
\item{breaks}{Content to be passed to the \code{breaks} argument of the generic \code{hist} function, if \code{histogram} is \code{TRUE}.  Default is 50, meaning that 50 vertical bars with equal binwidths will be drawn.} 
\item{col}{Colors of the plotting characters displaying the cluster assignment (if \code{labels} is \code{TRUE}).  If \code{NULL} (default), it will be determined automatically.}
\item{pch}{Plotting character used to show the cluster assignment.}
\item{cex}{Size of the plotting character showing the cluster assignment.}
\item{\dots}{Further arguments passed to \code{curve} (and also \code{hist} if \code{histogram} is \code{TRUE}).}
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{flowClust}}, \code{\link[=plot.flowClust]{plot}}, \code{\link[=density.flowClust]{density}}
}

\examples{
res <- flowClust(iris[,1:4], K=3)
hist(res, data=iris, subset="Petal.Length", breaks=30)
hist(res, data=iris, subset=3, histogram=FALSE, labels=FALSE)
hist(res, data=iris, subset="Petal.Length", breaks=30, include=2:3)
}

\keyword{graphs}





