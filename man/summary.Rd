\name{summary,flowClust-method}
\docType{methods}
\alias{summary,flowClust-method}
\alias{summary,flowClustList-method}
\alias{summary,tmixFilterResult-method}
\alias{summary,tmixFilterResultList-method}
\alias{summary.flowClust}
\alias{summary.tmixFilterResult}

\title{Summary Method for flowClust Object}

\description{
This method prints out various characteristics of the model fitted via robust model-based clustering.
}

\usage{
\S4method{summary}{flowClust}(object)

\S4method{summary}{tmixFilterResult}(object)
}

\arguments{
\item{object}{Object returned from \code{\link{flowClust}} or from \code{\link[=tmixFilter]{filter}}.}
}

\details{
Various characteristics of the fitted model will be given under the following five categories: Experiment Information, Clustering Summary, Transformation Parameter, Information Criteria, and Data Quality.  Under Data Quality, information about data filtering, outliers, and uncertainty is given.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\seealso{
\code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}}, \code{\link[=show,flowClust-method]{show}}
}

\keyword{print}
