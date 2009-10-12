\name{tmixFilter}
\docType{methods}
\alias{tmixFilter}
\alias{filter,flowFrame,tmixFilter-method}
\alias{filter,flowFrame,filter-method}
\alias{filter,flowFrame,filterSet-method}
\alias{filter,flowSet,filter-method}
\alias{filter,flowSet,filterSet-method}
\alias{filter,flowSet,list-method}
\alias{filter-method}
\alias{filter.flowFrame}
\alias{filter}
\alias{tmixFilter-class}
\alias{tmixFilterResult-class}
\alias{tmixFilterResultList}
\alias{tmixFilterResultList-class}
\alias{\%in\%,flowFrame,tmixFilter-method}
\alias{\%in\%,ANY,tmixFilterResult-method}
\alias{\%in\%,flowFrame,tmixFilterResult-method}
\alias{\%in\%,ANY,tmixFilterResultList-method}
\alias{[,flowFrame,tmixFilterResult-method}
\alias{[,flowFrame,tmixFilterResultList-method}
\alias{[[,tmixFilterResultList,ANY-method}

\title{Creating Filters and Filtering Flow Cytometry Data}

\description{
The \code{tmixFilter} function creates a filter object which is then passed to the \code{filter} method that performs filtering on a \code{flowFrame} object.  This method pair is provided to let \pkg{flowClust} integrate with the \pkg{flowCore} package.
}

\usage{
tmixFilter(filterId="tmixFilter", parameters="", \dots)

}

\arguments{
  
  \item{filterId}{A character string that identifies the filter created.}
  
  \item{parameters}{A character vector specifying the variables to be
    used in filtering.  When it is left unspecified, all the variables
    of the \code{flowFrame} object are used when running \code{filter}.
    Note that its content will be passed to the \code{varNames} argument
    of \code{\link{flowClust}} when running \code{filter}.}

  \item{\dots}{Other arguments passed to the \code{\link{flowClust}}
    function when running \code{\link[flowCore]{filter}},
    namely, \code{expName}, \code{K}, \code{B}, \code{tol}, \code{nu},
    \code{lambda}, \code{nu.est}, \code{trans}, \code{min.count}, \code{max.count},
    \code{min}, \code{max}, \code{level}, \code{u.cutoff},
    \code{z.cutoff}, \code{randomStart}, \code{B.init}, \code{tol.init},
    \code{seed} and \code{criterion}.  All arguments
    are optional except \code{K} that specifies the number of clusters.}
  
}

\value{
  
  The \code{tmixFilter} function returns an object of class
  \code{tmixFilter} that stores all the settings required for performing
  the \code{filter} operations.

  The \code{\link[flowCore]{filter}} method is defined in
  package \code{flowCore} and returns an object of class
  \code{tmixFilterResult} (or \code{tmixFilterResultList} if
  \code{filter@K} has a length >1) that stores the filtering results.  }


\note{
  
  The \code{tmixFilter} function returns an object of class
  \code{tmixFilter} that extends the virtual parent
  \code{\link[flowCore:filter-class]{filter}} class in the
  \pkg{flowCore} package.  Hence, the
  \link[flowCore:filter-class]{filter operators}, namely, \code{&},
  \code{|}, \code{!} and \code{\%subset\%}, also work for the
  \code{tmixFilter} class.


  If \code{filter@K} is of length 1, the \code{filter} method returns an
  object of class \code{tmixFilterResult}.  This class extends both the
  \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResult}}
  class (in the \pkg{flowCore} package) and the \code{\link{flowClust}}
  class.  Operations defined for the \code{multipleFilterResult} class,
  like \code{\link[flowCore:filter-class]{\%in\%}}, \code{\link{Subset}}
  and \code{\link{split}}, also work for the \code{tmixFilterResult}
  class.  Likewise, methods or functions designed to retrieve filtering
  (clustering) information from a \code{\link{flowClust}} object can
  also be applied on a \code{tmixFilterResult} object.  These include
  \code{\link{criterion}}, \code{\link{ruleOutliers}},
  \code{\link[=ruleOutliers]{ruleOutliers<-}}, \code{\link{Map}},
  \code{\link{posterior}}, \code{\link{importance}},
  \code{\link{uncertainty}} and \code{\link{getEstimates}}.  Various
  functionalities for plotting the filtering results are also available
  (see the links below).

  If \code{filter@K} has a length >1, the function returns an object of
  class \code{tmixFilterResultList}.  This class extends both the
  \code{\link{flowClustList}} class and the
  \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResult}}
  class.  Note that when a \code{tmixFilterResultList} object is used in
  place of a \code{tmixFilterResult} object, in most cases the list
  element corresponding to the best model will be extracted and passed
  to the method/function call.  }

\author{
Raphael Gottardo <\email{raph@stat.ubc.ca}>, Kenneth Lo <\email{c.lo@stat.ubc.ca}>
}

\references{
Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A} \bold{73}, 321-332.
}

\seealso{
\code{\link{flowClust}}, \code{\link[=summary.tmixFilterResult]{summary}}, \code{\link[=plot.tmixFilterResult]{plot}}, \code{\link[=density.flowClust]{density}}, \code{\link[=hist.flowClust]{hist}}, \code{\link{Subset}}, \code{\link{split}}, \code{\link{ruleOutliers}}, \code{\link{Map}}
}

\examples{
### The example below largely resembles the one in the flowClust
### man page.  The main purpose here is to demonstrate how the
### entire cluster analysis can be done in a fashion highly
### integrated into flowCore.


data(rituximab)

### create a filter object
s1filter <- tmixFilter("s1", c("FSC.H", "SSC.H"), K=1)
### cluster the data using FSC.H and SSC.H
res1 <- filter(rituximab, s1filter)

### remove outliers before proceeding to the second stage
# \%in\% operator returns a logical vector indicating whether each
# of the observations lies inside the gate or not
rituximab2 <- rituximab[rituximab \%in\% res1,]
# a shorthand for the above line
rituximab2 <- rituximab[res1,]
# this can also be done using the Subset method
rituximab2 <- Subset(rituximab, res1)

### cluster the data using FL1.H and FL3.H (with 3 clusters)
s2filter <- tmixFilter("s2", c("FL1.H", "FL3.H"), K=3)
res2 <- filter(rituximab2, s2filter)

show(s2filter)
show(res2)
summary(res2)

# to demonstrate the use of the split method
split(rituximab2, res2)
split(rituximab2, res2, population=list(sc1=c(1,2), sc2=3))

# to show the cluster assignment of observations
table(Map(res2))

# to show the cluster centres (i.e., the mean parameter estimates
# transformed back to the original scale) and proportions
getEstimates(res2)

### demonstrate the use of various plotting methods
# a scatterplot
plot(rituximab2, res2, level=0.8)
plot(rituximab2, res2, level=0.8, include=c(1,2), grayscale=TRUE,
    pch.outliers=2)
# a contour / image plot
res2.den <- density(res2, data=rituximab2)
plot(res2.den)
plot(res2.den, scale="sqrt", drawlabels=FALSE)
plot(res2.den, type="image", nlevels=100)
plot(density(res2, include=c(1,2), from=c(0,0), to=c(400,600)))
# a histogram (1-D density) plot
plot(rituximab2, res2, "FL1.H")

### to demonstrate the use of the ruleOutliers method
summary(res2)
# change the rule to call outliers
ruleOutliers(res2) <- list(level=0.95)
# augmented cluster boundaries lead to fewer outliers
summary(res2)

# the following line illustrates how to select a subset of data 
# to perform cluster analysis through the min and max arguments;
# also note the use of level to specify a rule to call outliers
# other than the default
s2t <- tmixFilter("s2t", c("FL1.H", "FL3.H"), K=3, B=100, 
    min=c(0,0), max=c(400,800), level=0.95, z.cutoff=0.5)
filter(rituximab2, s2t)
}

\keyword{cluster}
\keyword{models}
