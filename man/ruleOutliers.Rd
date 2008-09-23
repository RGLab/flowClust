\name{ruleOutliers,flowClust-method}
\docType{methods}
\alias{ruleOutliers,flowClust-method}
\alias{ruleOutliers,flowClustList-method}
\alias{ruleOutliers.flowClust}
\alias{ruleOutliers}
\alias{ruleOutliers<-,flowClust,list-method}
\alias{ruleOutliers<-,flowClustList,list-method}
\alias{ruleOutliers<-}

\title{Showing or Modifying the Rule used to Identify Outliers} 

\description{
This method shows or modifies the rule used to identify outliers.
}

\usage{
\S4method{ruleOutliers}{flowClust}(object)
ruleOutliers(object) <- value
}

\arguments{
\item{object}{Object returned from \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}}.}
\item{value}{A list object with one or more of the following named elements: \code{level}, \code{u.cutoff} and \code{z.cutoff}.  Their interpretations are the same as those of the corresponding arguments in the \code{\link{flowClust}} function.  Note that when both \code{level} and \code{u.cutoff} are missing, the rule set by the original value of \code{level} or \code{u.cutoff} will be unchanged rather than removed.  Likewise, when \code{z.cutoff} is missing, the rule set by the original value of \code{z.cutoff} will be retained.}
}

\value{
The replacement method modifies \code{object@ruleOutliers} (or \code{object[[k]]@ruleOutliers} if \code{object} is of class \code{flowClustList} or \code{tmixFilterResultList}) AND updates the logical vector \code{object@flagOutliers} (or \code{object[[k]]@ruleOutliers}) according to the new rule.
}

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}}
}

\examples{
res <- flowClust(iris[,1:4], K=3)
summary(res)
# change the rule to call outliers
ruleOutliers(res) <- list(z.cutoff=0.9)
# more outliers now
summary(res)
ruleOutliers(res) <- list(level=0.95, z.cutoff=0)
# augmented cluster boundaries lead to fewer outliers
summary(res)
}

\keyword{manip}
