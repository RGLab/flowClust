\name{flowClust}
\alias{flowClust}
\alias{flowClust-class}
\alias{flowClustList}
\alias{flowClustList-class}
\alias{\%in\%,ANY,flowClust-method}
\alias{\%in\%,ANY,flowClustList-method}
\alias{[,flowFrame,flowClust-method}
\alias{[,flowFrame,flowClustList-method}

\title{Robust Model-based Clustering for Flow Cytometry} 

\description{
This function performs automated clustering for identifying cell populations in flow cytometry data.  The approach is based on the \eqn{t} mixture model with the Box-Cox transformation, which provides a unified framework to handle outlier identification and data transformation simultaneously.
}

\usage{
flowClust(x, expName="Flow Experiment", varNames=NULL, K, B=500, 
          tol=1e-5, nu=4, lambda=1, nu.est=0, trans=1,
          min.count=10, max.count=10, min=NULL, max=NULL,
          level=0.9, u.cutoff=NULL, z.cutoff=0, randomStart=0, 
          B.init=B, tol.init=1e-2, seed=1, criterion="BIC",
          control=NULL)
}

\arguments{
\item{x}{A numeric vector, matrix, data frame of observations, or object of class \code{flowFrame}.  Rows correspond to observations and columns correspond to variables.}
\item{expName}{A character string giving the name of the experiment.}
\item{varNames}{A character vector specifying the variables (columns) to be included in clustering.  When it is left unspecified, all the variables will be used.}
\item{K}{An integer vector indicating the numbers of clusters.}
\item{B}{The maximum number of EM iterations.}
\item{tol}{The tolerance used to assess the convergence of the EM.}
\item{nu}{The degrees of freedom used for the \eqn{t} distribution.  Default is 4.  If \code{nu=Inf}, Gaussian distribution will be used.}
\item{lambda}{The initial transformation to be applied to the data.}
\item{nu.est}{A numeric indicating whether \code{nu} is to be estimated or not.  May take 0 (no estimation, default), 1 (estimation) or 2 (cluster-specific estimation).}
\item{trans}{A numeric indicating whether the Box-Cox transformation parameter is estimated from the data.  May take 0 (no estimation), 1 (estimation, default) or 2 (cluster-specific estimation).}
\item{min.count}{An integer specifying the threshold count for filtering data points from below.  The default is 10, meaning that if 10 or more data points are smaller than or equal to \code{min}, they will be excluded from the analysis.  If \code{min} is \code{NULL}, then the minimum of data as per each variable will be used.  To suppress filtering, set it as -1.}
\item{max.count}{An integer specifying the threshold count for filtering data points from above.  Interpretation is similar to that of \code{min.count}.}
\item{min}{The lower boundary set for data filtering.  Note that it is a vector of length equal to the number of variables (columns), implying that a different value can be set as per each variable.}
\item{max}{The upper boundary set for data filtering.  Interpretation is similar to that of \code{min}.}
\item{level}{A numeric value between 0 and 1 specifying the threshold quantile level used to call a point an outlier.  The default is 0.9, meaning that any point outside the 90\% quantile region will be called an outlier.}
\item{u.cutoff}{Another criterion used to identify outliers.  If this is \code{NULL}, then \code{level} will be used.  Otherwise, this specifies the threshold (e.g., 0.5) for \eqn{u}, a quantity used to measure the degree of \dQuote{outlyingness} based on the Mahalanobis distance.  Please refer to Lo et al. (2008) for more details.}
\item{z.cutoff}{A numeric value between 0 and 1 underlying a criterion which may be used together with \code{level}/\code{u.cutoff} to identify outliers.  A point with the probability of assignment \eqn{z} (i.e., the posterior probability that a data point belongs to the cluster assigned) smaller than \code{z.cutoff} will be called an outlier.  The default is 0, meaning that assignment will be made no matter how small the associated probability is, and outliers will be identified solely based on the rule set by \code{level} or \code{cutoff}.}
\item{randomStart}{A numeric value indicating how many times a random parition of the data is generated for initialization.  For instance, if \code{randomStart} is 10, 10 random partitions of the data will be generated, each of which is followed by a short EM run.  The partition leading to the highest likelihood value will be adopted to be the initial partition for the eventual long EM run.  The default is 0, meaning that this initialization strategy is not applied and hierarchical clustering is used instead.}
\item{B.init}{The maximum number of EM iterations following each random partition in random initialization.}
\item{tol.init}{The tolerance used as the stopping criterion for the short EM runs in random initialization.}
\item{seed}{An integer giving the seed number used when \code{randomStart>0}.}
\item{criterion}{A character string stating the criterion used to choose the best model.  May take either \code{"BIC"} or \code{"ICL"}.  This argument is only relevant when \code{length(K)>1}.}
\item{control}{An argument reserved for internal use.}
}

\details{
Estimation of the unknown parameters (including the Box-Cox parameter) is done via an Expectation-Maximization (EM) algorithm.  At each EM iteration, Brent's algorithm is used to find the optimal value of the Box-Cox transformation parameter.  Conditional on the transformation parameter, all other estimates can be obtained in closed form.  Please refer to Lo et al. (2008) for more details.

The \pkg{flowClust} package makes extensive use of the GSL as well as BLAS.  If an optimized BLAS library is provided when compiling the package, the \pkg{flowClust} package will be able to run multi-threaded processes.

Various operations have been defined for the object returned from \code{\link{flowClust}}.  These include:
\tabular{ll}{
Subsetting operations: \tab \code{\%in\%}, \code{\link{Subset}} and \code{\link{split}} \cr
Slot retrieval operations: \tab \code{\link{ruleOutliers}}, \code{\link{Map}}, \code{\link{criterion}}, \code{\link{posterior}}, \code{\link{importance}}, \code{\link{uncertainty}} and \code{\link{getEstimates}} \cr
Graphical operations: \tab \code{\link[=plot,flowClust-method]{plot}}, \code{\link[=density,flowClust-method]{density}} and \code{\link[=hist,flowClust-method]{hist}} \cr
}

In addition, to facilitate the integration with the \pkg{flowCore} package for processing flow cytometry data, the \code{flowClust} operation can be done through a method pair (\code{\link{tmixFilter}} and \code{\link[=tmixFilter]{filter}}) such that various methods defined in \pkg{flowCore} can be applied on the object created from the filtering operation.
}

\value{
If \code{K} is of length 1, the function returns an object of class \code{flowClust} containing the following slots, where \eqn{K} is the number of clusters, \eqn{N} is the number of observations and \eqn{P} is the number of variables:
\item{expName}{Content of the \code{expName} argument.}
\item{varNames}{Content of the \code{varNames} argument if provided; generated if available otherwise.}
\item{K}{An integer showing the number of clusters.}
\item{w}{A vector of length \eqn{K}, containing the estimates of the \eqn{K} cluster proportions.}
\item{mu}{A matrix of size \eqn{K \times P}{K x P}, containing the estimates of the \eqn{K} mean vectors.}
\item{sigma}{An array of dimension \eqn{K \times P \times P}{K x P x P}, containing the estimates of the \eqn{K} covariance matrices.}
\item{lambda}{The Box-Cox transformation parameter estimate.}
\item{nu}{The degrees of freedom for the \eqn{t} distribution.}
\item{z}{A matrix of size \eqn{N \times K}{N x K}, containing the posterior probabilities of cluster memberships.  The probabilities in each row sum up to one.}
\item{u}{A matrix of size \eqn{N \times K}{N x K}, containing the \dQuote{weights} (the contribution for computing cluster mean and covariance matrix) of each data point in each cluster.  Since this quantity decreases monotonically with the Mahalanobis distance, it can also be interpreted as the level of \dQuote{outlyingness} of a data point.  Note that, when \code{nu=Inf}, this slot is used to store the Mahalanobis distances instead.}
\item{label}{A vector of size \eqn{N}, showing the cluster membership according to the initial partition (i.e., hierarchical clustering if \code{randomStart} is \code{FALSE}).  Filtered observations will be labelled as \code{NA}.  Unassigned observations (which may occur since only 1500 observations at maximum are taken for hierarchical clustering) will be labelled as 0.}
\item{uncertainty}{A vector of size \eqn{N}, containing the uncertainty about the cluster assignment.  Uncertainty is defined as 1 minus the posterior probability that a data point belongs to the cluster to which it is assigned.}
\item{ruleOutliers}{A numeric vector of size 3, storing the rule used to call outliers.  The first element is 0 if the criterion is set by the \code{level} argument, or 1 if it is set by \code{u.cutoff}.  The second element copies the content of either the \code{level} or \code{u.cutoff} argument.  The third element copies the content of the \code{z.cutoff} argument.  For instance, if points are called outliers when they lie outside the 90\% quantile region or have assignment probabilities less than 0.5, then \code{ruleOutliers} is \code{c(0, 0.9, 0.5)}.  If points are called outliers only if their \dQuote{weights} in the assigned clusters are less than 0.5 regardless of the assignment probabilities, then \code{ruleOutliers} becomes \code{c(1, 0.5, 0)}.}
\item{flagOutliers}{A logical vector of size \eqn{N}, showing whether each data point is called an outlier or not based on the rule defined by \code{level}/\code{u.cutoff} and \code{z.cutoff}.}
\item{rm.min}{Number of points filtered from below.}
\item{rm.max}{Number of points filtered from above.}
\item{logLike}{The log-likelihood of the fitted mixture model.}
\item{BIC}{The Bayesian Information Criterion for the fitted mixture model.}
\item{ICL}{The Integrated Completed Likelihood for the fitted mixture model.}
If \code{K} has a length >1, the function returns an object of class \code{flowClustList}.  Its data part is a list with the same length as \code{K}, each element of which is a \code{flowClust} object corresponding to a specific number of clusters.  In addition, the resultant \code{flowClustList} object contains the following slots:\cr

\code{index} An integer giving the index of the list element corresponding to the best model as selected by \code{criterion}.\cr
\code{criterion} The criterion used to choose the best model -- either \code{"BIC"} or \code{"ICL"}.\cr

Note that when a \code{flowClustList} object is used in place of a \code{flowClust} object, in most cases the list element corresponding to the best model will be extracted and passed to the method/function call.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link[=summary.flowClust]{summary}}, \code{\link[=plot.flowClust]{plot}}, \code{\link[=density.flowClust]{density}}, \code{\link[=hist.flowClust]{hist}}, \code{\link{Subset}}, \code{\link{split}}, \code{\link{ruleOutliers}}, \code{\link{Map}}, \code{\link{SimulateMixture}}
}

\examples{
data(rituximab)

### cluster the data using FSC.H and SSC.H
res1 <- flowClust(rituximab, varNames=c("FSC.H", "SSC.H"), K=1)

### remove outliers before proceeding to the second stage
# \%in\% operator returns a logical vector indicating whether each
# of the observations lies within the cluster boundary or not
rituximab2 <- rituximab[rituximab \%in\% res1,]
# a shorthand for the above line
rituximab2 <- rituximab[res1,]
# this can also be done using the Subset method
rituximab2 <- Subset(rituximab, res1)

### cluster the data using FL1.H and FL3.H (with 3 clusters)
res2 <- flowClust(rituximab2, varNames=c("FL1.H", "FL3.H"), K=3)
show(res2)
summary(res2)

# to demonstrate the use of the split method
split(rituximab2, res2)
split(rituximab2, res2, population=list(sc1=c(1,2), sc2=3))

# to show the cluster assignment of observations
table(Map(res2))

# to show the cluster centres (i.e., the mean parameter estimates
# transformed back to the original scale)
getEstimates(res2)$locations

### demonstrate the use of various plotting methods
# a scatterplot
plot(res2, data=rituximab2, level=0.8)
# a contour / image plot
res2.den <- density(res2, data=rituximab2)
plot(res2.den)
plot(res2.den, type="image", nlevels=100)
# a histogram (1-D density) plot
hist(res2, data=rituximab2, subset="FL1.H")

# the following line illustrates how to select a subset of data 
# to perform cluster analysis through the min and max arguments;
# also note the use of level to specify a rule to call outliers
# other than the default
flowClust(rituximab2, varNames=c("FL1.H", "FL3.H"), K=3, B=100, 
    min=c(0,0), max=c(400,800), level=0.95, z.cutoff=0.5)
}

\keyword{cluster}
\keyword{models}
