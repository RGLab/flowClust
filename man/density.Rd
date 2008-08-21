\name{density,flowClust-method}
\docType{methods}
\alias{density-method}
\alias{density,flowClust-method}
\alias{flowDens-class}
\alias{density.flowClust}

\title{Grid of Density Values for the Fitted t Mixture Model with Box-Cox Transformation} 

\description{
This method constructs the \code{flowDens} object which is used to generate a contour or image plot.
}

\usage{
\S4method{density}{flowClust}(x, data=NULL, subset=c(1,2), include=1:(x@K), 
        npoints=c(100,100), from=NULL, to=NULL)
}

\arguments{
\item{x}{Object returned from \code{\link{flowClust}} or from running \code{\link[=filter.flowFrame]{filter}} on a \code{flowFrame} object.}
\item{data}{A matrix, data frame of observations, or object of class \code{flowFrame}.  This is the object on which \code{flowClust} or \code{filter} was performed.  If this argument is not specified, the grid square upon which densities will be computed must be provided (through arguments \code{from} and \code{to}).}
\item{subset}{A numeric vector of length two indicating which two variables are selected for the scatterplot.  Alternatively, a character vector containing the names of the two variables is allowed if \code{x@varNames} is not \code{NULL}.}
\item{include}{A numeric vector specifying which clusters are included to compute the density values.  By default, all clusters are included.}
\item{npoints}{A numeric vector of size two specifying the number of grid points in \eqn{x} (horizontal) and \eqn{y} (vertical) directions respectively.}
\item{from}{A numeric vector of size two specifying the coordinates of the lower left point of the grid square.  Note that, if this (and \code{to}) is not specified, \code{data} must be provided such that the range in the two variables (dimensions) selected will be used to define the grid square.}
\item{to}{A numeric vector of size two specifying the co-ordinates of the upper right point of the grid square.}
}

\details{
The \code{flowDens} object returned is to be passed to the \code{plot} method for generating a contour or image plot.
}

\value{
An object of class \code{flowDens} containing the following slots is constructed:
\item{dx}{A numeric vector of length \code{npoints[1]}; the \eqn{x}-coordinates of the grid points.}
\item{dy}{A numeric vector of length \code{npoints[2]}; the \eqn{y}-coordinates of the grid points.}
\item{value}{A matrix of size \code{npoints[1]} \eqn{\times}{x} \code{npoints[2]}; the density values at the grid points.}
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\seealso{
\code{\link[=plot.flowDens]{plot}}, \code{\link{flowClust}}
}

\examples{
res <- flowClust(iris[,1:4], K=3)

res.den <- density(res, data=iris)
plot(res.den)
plot(res.den, scale="sqrt", drawlabels=FALSE)
plot(res.den, type="image", nlevels=100)
plot(res.den, type="image", scale="sqrt", nlevels=100,
    color="gray")

plot(density(res, data=iris, subset=c("Petal.Length",
    "Petal.Width"), include=c(2,3)))
plot(density(res, subset=c(3,4), include=c(2,3), from=c(3,0.75),
    to=c(7,2.5)))
}

\keyword{graphs}





