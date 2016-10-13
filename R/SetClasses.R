#' @rdname flowClust
#' @importClassesFrom Biobase AssayData
#' @export 
setClass("flowClust",
         representation(expName="character", varNames="character",
                        K="numeric", w="vector", mu="matrix",
                        sigma="array", lambda="numeric", nu="numeric",
                        z="matrix", u="matrix", label="vector",
                        uncertainty="vector", ruleOutliers="vector",
                        flagOutliers="vector", rm.min="numeric",
                        rm.max="numeric", logLike="numeric", BIC="numeric",
                        ICL="numeric", usePrior="character", prior="list"),
         prototype(expName=character(0), varNames=character(0), K=numeric(0),
                   w=rep(numeric(0),0), mu=matrix(numeric(0), nrow=0, ncol=0),
                   sigma=array(numeric(0), c(0,0,0)), lambda=numeric(0),
                   nu=numeric(0), z=matrix(numeric(0), nrow=0, ncol=0),
                   u=matrix(numeric(0), nrow=0, ncol=0),
                   label=rep(numeric(0),0), uncertainty=rep(numeric(0),0),
                   ruleOutliers=rep(numeric(0),0),
                   flagOutliers=rep(logical(0), 0), rm.min=numeric(0),
                   rm.max=numeric(0), logLike=numeric(0), BIC=numeric(0),
                   ICL=numeric(0), usePrior="no", prior=list(NA)))


#' @rdname flowClust
#' @export 
setClass("flowClustList",
         representation("list", index="numeric", criterion="character"),
         prototype(vector("list",0), index=numeric(0), criterion="BIC"),
         validity=function(object)
                  {
                      if (!all(sapply(object@.Data, is, class2="flowClust")))
                          return("Not a list of flowClust results!")
                      else return(TRUE)
                  })

#' @rdname density
#' @export 
setClass("flowDens",
         representation(dx="matrix", dy="matrix", value="matrix"),
         prototype(dx=matrix(numeric(0), nrow=0, ncol=1),
                   dy=matrix(numeric(0), nrow=0, ncol=1),
                   value=matrix(numeric(0), nrow=0, ncol=0)))

#' @rdname tmixFilter
#' @export 
setClass("tmixFilter",
         representation(expName="character", K="numeric", B="numeric",
                        tol="numeric", nu="numeric", lambda="numeric",
                        nu.est="numeric", trans="numeric", min.count="numeric",
                        max.count="numeric", min="vector", max="vector",
                        level="numeric", u.cutoff="numeric", z.cutoff="numeric",
                        randomStart="numeric", B.init="numeric", tol.init="numeric",
                        seed="numeric", criterion="character", control="list",
                        usePrior="character", prior="list"),
         prototype(expName="Flow Experiment", K=numeric(0), B=500, tol=1e-5,
                   nu=4, lambda=1, nu.est=0, trans=1, min.count=10, max.count=10,
                   min=NA, max=NA, level=0.9, u.cutoff=NA_real_, z.cutoff=0,
                   randomStart=0, B.init=500, tol.init=1e-2, seed=1, criterion="BIC", 
                   control=vector("list",0), usePrior="no", prior=list(NA)),
         contains="parameterFilter")




#' Creating Filters and Filtering Flow Cytometry Data
#' 
#' The \code{tmixFilter} function creates a filter object which is then passed
#' to the \code{filter} method that performs filtering on a \code{flowFrame}
#' object.  This method pair is provided to let \pkg{flowClust} integrate with
#' the \pkg{flowCore} package.
#' 
#' @name tmixFilter
#' 
#' @param filterId A character string that identifies the filter created.
#' @param parameters A character vector specifying the variables to be used in
#' filtering.  When it is left unspecified, all the variables of the
#' \code{flowFrame} object are used when running \code{filter}.  Note that its
#' content will be passed to the \code{varNames} argument of
#' \code{\link{flowClust}} when running \code{filter}.
#' 
#' @param \dots Other arguments passed to the \code{\link{flowClust}} function
#' when running \code{\link[flowCore]{filter}}, namely, \code{expName},
#' \code{K}, \code{B}, \code{tol}, \code{nu}, \code{lambda}, \code{nu.est},
#' \code{trans}, \code{min.count}, \code{max.count}, \code{min}, \code{max},
#' \code{level}, \code{u.cutoff}, \code{z.cutoff}, \code{randomStart},
#' \code{B.init}, \code{tol.init}, \code{seed} and \code{criterion}.  All
#' arguments are optional except \code{K} that specifies the number of
#' clusters.
#' 
#' The \code{tmixFilter} function returns an object of class \code{tmixFilter}
#' that stores all the settings required for performing the \code{filter}
#' operations.
#' 
#' The \code{\link[flowCore]{filter}} method is defined in package
#' \code{flowCore} and returns an object of class \code{tmixFilterResult} (or
#' \code{tmixFilterResultList} if \code{filter} has a length >1) that stores
#' the filtering results.
#'  
#' The \code{tmixFilter} function returns an object of class \code{tmixFilter}
#' that extends the virtual parent \code{\link[flowCore:filter-class]{filter}}
#' class in the \pkg{flowCore} package.  Hence, the
#' \link[flowCore:filter-class]{filter operators}, namely, \code{&}, \code{|},
#' \code{!} and \code{subset}, also work for the \code{tmixFilter} class.
#' 
#' If \code{filter} is of length 1, the \code{filter} method returns an
#' object of class \code{tmixFilterResult}.  This class extends both the
#' \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResult}}
#' class (in the \pkg{flowCore} package) and the \code{\link{flowClust}} class.
#' Operations defined for the \code{multipleFilterResult} class, like
#' \code{\link[flowCore:filter-class]{\%in\%}}, \code{\link{Subset}} and
#' \code{\link{split}}, also work for the \code{tmixFilterResult} class.
#' Likewise, methods or functions designed to retrieve filtering (clustering)
#' information from a \code{\link{flowClust}} object can also be applied on a
#' \code{tmixFilterResult} object.  These include \code{\link{criterion}},
#' \code{\link{ruleOutliers}}, \code{\link[=ruleOutliers]{ruleOutliers<-}},
#' \code{\link{Map}}, \code{\link{posterior}}, \code{\link{importance}},
#' \code{\link{uncertainty}} and \code{\link{getEstimates}}.  Various
#' functionalities for plotting the filtering results are also available (see
#' the links below).
#' 
#' If \code{filter} has a length >1, the function returns an object of class
#' \code{tmixFilterResultList}.  This class extends both the
#' \code{\link{flowClustList}} class and the
#' \code{\link[flowCore:multipleFilterResult-class]{multipleFilterResult}}
#' class.  Note that when a \code{tmixFilterResultList} object is used in place
#' of a \code{tmixFilterResult} object, in most cases the list element
#' corresponding to the best model will be extracted and passed to the
#' method/function call.
#' 
#' @seealso \code{\link{flowClust}},
#' \code{\link[=summary.tmixFilterResult]{summary}},
#' \code{\link[=plot.tmixFilterResult]{plot}},
#' \code{\link[=density.flowClust]{density}},
#' \code{\link[=hist.flowClust]{hist}}, \code{\link{Subset}},
#' \code{\link{split}}, \code{\link{ruleOutliers}}, \code{\link{Map}}
#' 
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords cluster models
#' @examples
#' 
#' ### The example below largely resembles the one in the flowClust
#' ### man page.  The main purpose here is to demonstrate how the
#' ### entire cluster analysis can be done in a fashion highly
#' ### integrated into flowCore.
#' 
#' 
#' data(rituximab)
#' 
#' ### create a filter object
#' s1filter <- tmixFilter("s1", c("FSC.H", "SSC.H"), K=1)
#' ### cluster the data using FSC.H and SSC.H
#' res1 <- filter(rituximab, s1filter)
#' 
#' ### remove outliers before proceeding to the second stage
#' # %in% operator returns a logical vector indicating whether each
#' # of the observations lies inside the gate or not
#' rituximab2 <- rituximab[rituximab %in% res1,]
#' # a shorthand for the above line
#' rituximab2 <- rituximab[res1,]
#' # this can also be done using the Subset method
#' rituximab2 <- Subset(rituximab, res1)
#' 
#' ### cluster the data using FL1.H and FL3.H (with 3 clusters)
#' s2filter <- tmixFilter("s2", c("FL1.H", "FL3.H"), K=3)
#' res2 <- filter(rituximab2, s2filter)
#' 
#' show(s2filter)
#' show(res2)
#' summary(res2)
#' 
#' # to demonstrate the use of the split method
#' split(rituximab2, res2)
#' split(rituximab2, res2, population=list(sc1=c(1,2), sc2=3))
#' 
#' # to show the cluster assignment of observations
#' table(Map(res2))
#' 
#' # to show the cluster centres (i.e., the mean parameter estimates
#' # transformed back to the original scale) and proportions
#' getEstimates(res2)
#' 
#' ### demonstrate the use of various plotting methods
#' # a scatterplot
#' plot(rituximab2, res2, level=0.8)
#' plot(rituximab2, res2, level=0.8, include=c(1,2), grayscale=TRUE,
#'     pch.outliers=2)
#' # a contour / image plot
#' res2.den <- density(res2, data=rituximab2)
#' plot(res2.den)
#' plot(res2.den, scale="sqrt", drawlabels=FALSE)
#' plot(res2.den, type="image", nlevels=100)
#' plot(density(res2, include=c(1,2), from=c(0,0), to=c(400,600)))
#' # a histogram (1-D density) plot
#' plot(rituximab2, res2, "FL1.H")
#' 
#' ### to demonstrate the use of the ruleOutliers method
#' summary(res2)
#' # change the rule to call outliers
#' ruleOutliers(res2) <- list(level=0.95)
#' # augmented cluster boundaries lead to fewer outliers
#' summary(res2)
#' 
#' # the following line illustrates how to select a subset of data 
#' # to perform cluster analysis through the min and max arguments;
#' # also note the use of level to specify a rule to call outliers
#' # other than the default
#' s2t <- tmixFilter("s2t", c("FL1.H", "FL3.H"), K=3, B=100, 
#'     min=c(0,0), max=c(400,800), level=0.95, z.cutoff=0.5)
#' filter(rituximab2, s2t)
#' 
#' @rdname tmixFilter
#' @export 
tmixFilter <- function(filterId="tmixFilter", parameters="", ...)
{
    if (!"K" %in% names(list(...)))
        stop("The number of clusters, K, must be provided.")
    new("tmixFilter", filterId=filterId, parameters=as.list(parameters), ...)
}


#' @export 
#' @rdname tmixFilter
setClass("tmixFilterResult", contains=c("flowClust", "multipleFilterResult"))

#' @name tmixFilterResultList
#' @exportClass  tmixFilterResultList
#' @rdname tmixFilter
setClassUnion(name = "tmixFilterResultList",
    members=c("flowClustList", "multipleFilterResult"))
