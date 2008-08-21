\name{plot,flowClust-method}
\docType{methods}
\alias{plot,flowClust,missing-method}
\alias{plot,flowClust-method}
\alias{plot.flowClust}

\title{Scatterplot of Clustering Results} 

\description{
This method generates scatterplot revealing the cluster assignment, cluster boundaries according to the specified percentile as well as supplemental information like outliers or filtered observations.
}

\usage{
\S4method{plot}{flowClust}(x, data, subset=c(1,2), ellipse=TRUE, show.outliers=TRUE,
     show.rm=FALSE, include=1:(x@K), main=NULL, grayscale=FALSE,
     col=(if (grayscale) gray(1/4) else 2:(length(include)+1)),
     pch=".", cex=0.6, col.outliers=gray(3/4), pch.outliers=".",
     cex.outliers=cex, col.rm=1, pch.rm=1, cex.rm=0.6, ecol=1, 
     elty=1, level=NULL, u.cutoff=NULL, z.cutoff=NULL, 
     npoints=501, add=FALSE, \dots)
}

\arguments{
\item{x}{Object returned from \code{\link{flowClust}}.}
\item{data}{A matrix, data frame of observations, or object of class \code{flowFrame}. This is the object on which \code{flowClust} was performed.}
\item{subset}{A numeric vector of length two indicating which two variables are selected for the scatterplot.  Alternatively, a character vector containing the names of the two variables is allowed if \code{x@varNames} is not \code{NULL}.}
\item{ellipse}{A logical value indicating whether the cluster boundary is to be drawn or not.  If \code{TRUE}, the boundary will be drawn according to the level specified by \code{level} or \code{cutoff}.}
\item{show.outliers}{A logical value indicating whether outliers will be explicitly shown or not.}
\item{show.rm}{A logical value indicating whether filtered observations will be shown or not.}
\item{include}{A numeric vector specifying which clusters will be shown on the plot.  By default, all clusters are included.}
\item{main}{Title of the plot.}
\item{grayscale}{A logical value specifying if a grayscale plot is desired.  This argument takes effect only if the default values of relevant graphical arguments are taken.}
\item{col}{Color(s) of the plotting characters.  May specify a different color for each cluster.}
\item{pch}{Plotting character(s) of the plotting characters.  May specify a different character for each cluster.} 
\item{cex}{Size of the plotting characters.  May specify a different size for each cluster.}
\item{col.outliers}{Color of the plotting characters denoting outliers.}
\item{pch.outliers}{Plotting character(s) used to denote outliers.  May specify a different character for each cluster.}
\item{cex.outliers}{Size of the plotting characters used to denote outliers.  May specify a different size for each cluster.}
\item{col.rm}{Color of the plotting characters denoting filtered observations.}
\item{pch.rm}{Plotting character used to denote filtered observations.}
\item{cex.rm}{Size of the plotting character used to denote filtered observations.}
\item{ecol}{Color(s) of the lines representing the cluster boundaries.  May specify a different color for each cluster.}
\item{elty}{Line type(s) drawing the cluster boundaries.  May specify a different line type for each cluster.}
\item{level, u.cutoff, z.cutoff}{These three optional arguments specify the rule used to identify outliers.  By default, all of them are left unspecified, meaning that the rule stated in \code{x@ruleOutliers} will be taken.  Otherwise, these arguments will be passed to \code{\link{ruleOutliers}}.}
\item{npoints}{The number of points used to draw each cluster boundary.}
\item{add}{A logical value.  If \code{TRUE}, add to the current plot.}
\item{\dots}{Further graphical parameters passed to the generic function \code{plot}.}
}

\note{
The cluster boundaries need not be elliptical since Box-Cox transformation has been performed.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{flowClust}}
}

\examples{
res <- flowClust(iris[,1:4], K=3)

plot(res, data=iris, subset=c("Sepal.Width", "Petal.Width"),
    pch=1, pch.outlier=2)
plot(res, data=iris, subset=c(2,4), ellipse=FALSE,
    show.outlier=FALSE, pch=1)

plot(res, data=iris, subset=c("Petal.Length", "Petal.Width"),
    include=c(2,3), pch=1, pch.outlier=1)
plot(res, data=iris, subset=c("Petal.Length", "Petal.Width"),
    include=c(2,3), grayscale=TRUE, pch=1, pch.outlier=1,
    level=0.8)
}

\keyword{graphs}





