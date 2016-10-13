\name{flowClust-package}
\alias{flowClust-package}
\docType{package}
\title{
Clustering for Flow Cytometry
}
\description{
Robust model-based clustering using a \eqn{t} mixture model with Box-Cox transformation.
}
\details{
\tabular{ll}{
Package: \tab flowClust\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Depends: \tab R(>= 2.5.0), methods, mnormt, mclust, ellipse, Biobase, flowCore\cr
Collate: \tab SetClasses.R SetMethods.R plot.R flowClust.R SimulateMixture.R\cr
biocViews: \tab Clustering, Statistics, Visualization\cr
License: \tab Artistic-2.0\cr
Built: \tab R 2.6.1; universal-apple-darwin8.10.1; 2008-03-26 20:54:42; unix\cr
}
}

\section{Index}{
\describe{
\item{\code{\link{box}}}{Box-Cox Transformation}
\item{\code{\link{density,flowClust-method}}}{Grid of Density Values for the Fitted \eqn{t} Mixture Model with Box-Cox Transformation}
\item{\code{\link{dmvt}}}{Density of the Multivariate \eqn{t} Distribution with Box-Cox Tranformation}
\item{\code{\link{dmvtmix}}}{Density of the Multivariate \eqn{t} Mixture Distribution with Box-Cox Tranformation}
\item{\code{\link{flowClust}}}{Robust Model-based Clustering for Flow Cytometry}
\item{\code{\link{hist.flowClust}}}{1-D Density Plot (Histogram) of Clustering Results}
\item{\code{\link{Map,flowClust-method}}}{Cluster Assignment Based on Clustering Results}
\item{\code{\link[=posterior]{miscellaneous}}}{Various Functions for Retrieving Information from Clustering Results}
\item{\code{\link{plot,flowClust-method}}}{Scatterplot of Clustering Results}
\item{\code{\link{plot,flowDens-method}}}{Contour or Image Plot of Clustering Results}
\item{\code{\link{plot,flowFrame,tmixFilterResult-method}}}{Scatterplot / 1-D Density Plot of Filtering (Clustering) Results}
\item{\code{\link{rbox}}}{Reverse Box-Cox Transformation}
\item{\code{\link{ruleOutliers,flowClust-method}}}{Showing or Modifying the Rule used to Identify Outliers}
\item{\code{\link{show,flowClust-method}}}{Show Method for \code{flowClust} / \code{tmixFilterResult} Object}
\item{\code{\link{show,tmixFilter-method}}}{Show Method for \code{tmixFilter} Object}
\item{\code{\link{SimulateMixture}}}{Random Generation from a \eqn{t} Mixture Model with Box-Cox Transformation}
\item{\code{\link{split,flowClust-method}}}{Splitting Data Based on Clustering Results}
\item{\code{\link{Subset,flowClust-method}}}{Subsetting Data Based on Clustering Results}
\item{\code{\link{summary,flowClust-method}}}{Summary Method for \code{flowClust} Object}
\item{\code{\link{tmixFilter}}}{Creating Filters and Filtering Flow Cytometry Data}
}
}

\note{Further information is available in the vignette.}

\author{
Raphael Gottardo <raph@stat.ubc.ca>, Kenneth Lo <c.lo@stat.ubc.ca>

Maintainer: Raphael Gottardo <raph@stat.ubc.ca>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\keyword{package}
