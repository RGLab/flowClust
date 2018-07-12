#' @include flowClust.R
NULL

#' 1-D Density Plot (Histogram) of Clustering Results
#' 
#' This method generates a one-dimensional density plot for the specified
#' dimension (variable) based on the robust model-based clustering results.  A
#' histogram of the actual data or cluster assignment is optional for display.
#' 
#' @method hist flowClust
#' @param x Object returned from \code{\link{flowClust}} or from running
#' \code{filter} on a \code{flowFrame} object.
#' @param data A numeric vector, matrix, data frame of observations, or object
#' of class \code{flowFrame}. This is the object on which \code{flowClust} or
#' \code{filter} was performed.
#' @param subset An integer indicating which variable is selected for the plot.
#' Alternatively, a character string containing the name of the variable is
#' allowed if \code{x@varNames} is not \code{NULL}.
#' @param include A numeric vector specifying which clusters are shown on the
#' plot.  By default, all clusters are included.
#' @param histogram A logical value indicating whether a histogram of the
#' actual data is made in addition to the density plot or not.
#' @param labels A logical value indicating whether information about cluster
#' assignment is shown or not.
#' @param ylab Labels for the \eqn{x}- and \eqn{y}-axes respectively.
#' @param main Title of the plot.
#' @param col Colors of the plotting characters displaying the cluster
#' assignment (if \code{labels} is \code{TRUE}).  If \code{NULL} (default), it
#' will be determined automatically.
#' @param pch Plotting character used to show the cluster assignment.
#' @param cex Size of the plotting character showing the cluster assignment.
#' @param \dots other arguments
#' 
#'    xlim The range of \eqn{x}-values for the plot.  If \code{NULL}, the data range will be used.
#' 
#'    ylim The range of \eqn{y}-values for the plot.  If \code{NULL}, an optimal range will be determined automatically.
#' 
#'    breaks Content to be passed to the \code{breaks} argument of the generic \code{hist} function, if \code{histogram} is \code{TRUE}.  Default
#'            is 50, meaning that 50 vertical bars with equal binwidths will be drawn.
#' 
#'    ... Further arguments passed to \code{curve} (and also \code{hist}
#' 
#' if \code{histogram} is \code{TRUE}).
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}, \code{\link[=plot,flowClust-method]{plot}},
#' \code{\link[=density.flowClust]{density}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords graphs
#' @rdname hist
#' @export 
hist.flowClust <- function(x, data=NULL, subset=1, include=1:(x@K)
                           , histogram=TRUE
                           , labels=TRUE
                           , ylab="Density"
                           , main=NULL
                           , col=NULL, pch=20, cex=0.6
                           , ...)
{
  
  pre.obj <- .hist.flowClust(x = x, data = data, subset = subset, include = include, histogram = histogram, labels = labels, ...)

  den <- function(y) {
    flowClust.den(x = y, obj = x, subset = pre.obj$subset, include = include)
  }
  
  if (histogram) 
    hist(pre.obj$data2, breaks=pre.obj$tbreaks, freq=F
         , xlim=pre.obj$xlim, ylim=pre.obj$ylim, xlab=pre.obj$xlab, ylab=ylab, main=main
         , ...)
  
  curve(den, add=histogram
        , xlim=pre.obj$xlim
        , ylim=pre.obj$ylim
        , xlab=pre.obj$xlab
        , ylab=ylab, main=main
        , ...)
  
  if (labels) {
    if (is.null(col)) {
      if (length(include)<=4) col <- c("red", "blue", "green", "black")  else col <- 2:(length(include)+1)
    } else col<-matrix(col, length(include))
    j <- 0
    for (k in include) stripchart(pre.obj$data[Map(x, rm.outliers=F)==k], add=T
                                  , at=pre.obj$ymin - (pre.obj$ylim[2]-pre.obj$ymin)/100*(j<-j+1)
                                  , pch=pch, cex=cex, col=col[j])
  }
}

#' preprocessing flowClust results to prepare for the hist plot
#' It is helpful to separate this logic from the default hist plot function
#' so that it can be reused by the other kind of plot engine (e.g. ggplot)
#' @noRd
.hist.flowClust <- function(x, data=NULL, subset=1, include=1:(x@K)
                            , histogram=TRUE
                            , labels=TRUE, xlim=NULL, ylim=NULL
                            , xlab=(if (is.numeric(subset)) NULL else subset)
                            , breaks=50
                            , ...)
{
  
  
  if (is(data, "flowFrame")) data <- exprs(data)[,x@varNames,drop=FALSE]  else
    if (is(data, "matrix")) (if (length(x@varNames)>0) data <- as.matrix(data[,x@varNames]))  else
      if (is(data, "data.frame")) data <- as.matrix(data[,x@varNames])  else
        if (is(data, "vector")) data <- matrix(data)
        
        if (is.null(xlab) && x@varNames!="Not Available") xlab <- x@varNames[subset]
        if (!is.numeric(subset)) subset <- match(subset, x@varNames)
        # look for highest density value
        data <- data[,subset]
        data1 <- data[!is.na(x@flagOutliers)]
        if (is.null(ylim)) {
          tseq <- seq(min(data1), max(data1), length.out=500)
          ymax <- max(flowClust.den(x= tseq[tseq!=0], obj = x, subset = subset, include = include))
        }
        
        # look for highest point in histogram
        if (histogram) {
          data2 <- data[!is.na(x@flagOutliers) & is.element(Map(x, rm.outliers=F), include)]
          tbreaks <- hist(data1, breaks=breaks, plot=F)$breaks
          if (is.null(ylim)) {
            tplot <- hist(data2, breaks=tbreaks, plot=F)
            ymax <- max(ymax, tplot$density)
          }
        }
        
        if (is.null(xlim)) xlim <- range(data1)
        if (is.null(ylim)) ylim <- c(0, ymax)
        ymin <- ylim[1]
        if (labels) ylim[1] <- ylim[1] - (ylim[2]-ylim[1])/100*length(include)
        
        list(data = data, data2 = data2, xlim = xlim, ylim = ylim, ymax = ymax, ymin = ymin, xlab = xlab, tbreaks = tbreaks, subset = subset)
}
#' generate the curve that reflects the tmixture fitting outcome
#' 
#' @param x the numeric vector represents the x coordinates in plot
#' @param obj the flowClust object
#' @inheritParams hist.flowClust
flowClust.den <- function(x, obj, subset, include) {
  value <- 0
  nu <- rep(obj@nu, length.out=obj@K)
  
  if (length(obj@lambda)>0&(any(obj@lambda!=1))) {
    lambda <- rep(obj@lambda, length.out=obj@K)
    for (k in include) {
      xTrans <- (sign(x)*abs(x)^lambda[k] - 1) / lambda[k]
      sigma <- obj@sigma[k,subset,subset]
      if(!is.na(sigma))
        value <- value + obj@w[k] * dmvt(xTrans
                                         , obj@mu[k,subset]
                                         , sigma
                                         , nu[k]
                                         , log=F)$value * abs(x)^(lambda[k]-1)
    }
  } else {
    for (k in include){
      sigma <- obj@sigma[k,subset,subset]
      if(!is.na(sigma))
        value <- value + obj@w[k] * dmvt(x, obj@mu[k,subset], sigma, nu[k], log=F)$value
    } 
    
  }
  value <- value / sum(obj@w[include])
  value
}

#' @method hist flowClustList
#' @rdname hist
#' @export 
hist.flowClustList <- function(x, ...)
{
  x <- as(x, "flowClust")
  hist(x, ...)  
}

