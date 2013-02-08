\name{plotPrior}
\alias{plotPrior}
\title{
Plots a flowClust prior over some data.
}
\description{
Plots a flowClust prior overlaid on data.
}
\usage{
plotPrior(data, prior, dims = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
On object of class "flowFrame". The data to be plotted.
}
  \item{prior}{
An object of class "flowClustPrior", or "flowClustPriorList", returned by a call to \code{mkPrior}.
}
  \item{dims}{
A character vector of the dimensions to be included in  the plot. The dimension names should match column names in the prior and in the flowFrame.
}
  \item{\dots}{
Additional arguments to plotting functions, such as \code{smooth=TRUE/FALSE}

}
}
\details{
Generates a plot of a "flowClustPrior" or "flowClustPriorList" object overlaid on some data. Plots the prior means (Mu0), prior covariance of the means (Omega0), and prior sample covariance (Lambda0).
}
\value{
Silently returns zero.
}
\author{
Greg Finak <gfinak@fhcrc.org>
}


\keyword{ aplot }
\keyword{ dplot }% __ONLY ONE__ keyword per line
