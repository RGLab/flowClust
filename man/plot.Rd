\name{plot}
\alias{plot}
\alias{plot-methods}
\alias{plot,flowClustPriorTree,GatingHierarchy-method}
\alias{plot,flowClustPriorTree,GatingSet-method}
\title{
Plot a prior from a flowClustPriorTree
}
\description{
Plot a prior from a flowClustPriorTree.
}
\usage{
}
\arguments{
  \item{x}{
	\code{flowClustPriorTree}
}
  \item{y}{
	\code{GatingSet} or \code{GatingHierarchy}
}
\item{node}{
\code{numeric} The index of the node to be plotted.
}
\item{parent}{
\code{logical} if TRUE, plot the data from the parent. if FALSE, plot just the data that generated the prior.
}
\code{jitter}{
\code{logical} Add jitter. Useful for visualizaiton in some cases
}
\item{add}{
\code{logical} TRUE or FALSE, indicating whether to add the priors to the current plot.
}
\item{dims}{
Default \code{NULL}: dimensions are taken from the prior. \code{character} vector of the two dimensions to be plotted.
}
  \item{\dots}{
	Additional arguments to plot().
}
}
\details{
Plot a prior distribution, including prior mean, prior covariance, and covariance of the means. If \code{y} is a GatingHierarchy, the data from that sample will be plotted. If \code{y} is a GatingSet, the data from all samples in the GatingSet will be plotted, as well as the means and covariances of each sample. \code{parent=TRUE|FALSE} controls whether the data from the parent population is plotted.
}
\references{
\url{http://www.rglab.org/}
}

\author{
Greg Finak \email{gfinak@fhcrc.org}
}


\examples{
\dontrun{
#G is a GatingHierarchy
#g is a flowClustPriorTree
plot(g,G,3) #plot the third node
}
}
\keyword{ hplot }
