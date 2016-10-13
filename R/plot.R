#' Box-Cox Transformation
#' 
#' This function performs Box-Cox transformation on the inputted data matrix.
#' 
#' To allow for negative data values, a slightly modified version of the
#' original Box-Cox (1964) is used here.  This modified version originated from
#' Bickel and Doksum (1981), taking the following form: \deqn{f(y) =
#' \frac{\mathrm{sgn}(y)|y|^\lambda-1}{\lambda}}{f(y) = ( sgn(y)
#' abs(y)^(lambda) -1 ) / lambda} When negative data values are involved, the
#' transformation parameter, \eqn{\lambda}{\code{lambda}}, has to be positive
#' in order to avoid discontinuity across zero.
#' 
#' @param data A numeric vector, matrix or data frame of observations.
#' Negative data values are permitted.
#' @param lambda The transformation to be applied to the data.  If negative
#' data values are present, \code{lambda} has to be positive.
#' @return A numeric vector, matrix or data frame of the same dimension as
#' \code{data} is returned.
#' @seealso \code{\link{rbox}}
#' @references Bickel, P. J. and Doksum, K. A. (1981) An Analysis of
#' Transformations Revisited. \emph{J. Amer. Statist. Assoc.} \bold{76}(374),
#' 296-311.
#' 
#' Box, G. E. P. and Cox, D. R. (1964) An Analysis of Transformations. \emph{J.
#' R. Statist. Soc. B} \bold{26}, 211-252.
#' @keywords math
#' @examples
#' 
#' data(rituximab)
#' data <- exprs(rituximab)
#' summary(data)
#' # Transform data using Box-Cox with lambda=0.3
#' dataTrans <- box(data, 0.3)
#' # Reverse transform data; this should return back to the original rituximab data
#' summary(rbox(dataTrans, 0.3))
#' @export 
box <- function(data, lambda) {
    if (length(lambda)>1 || lambda!=0) data <- (sign(data)*abs(data)^lambda-1)/lambda else data <- log(data)
    data
}

#' Reverse Box-Cox Transformation
#' 
#' This function performs back transformation on Box-Cox transformed data.
#' 
#' 
#' @param data A numeric vector, matrix or data frame of observations.
#' @param lambda The Box-Cox transformation applied which results in the
#' inputted data matrix.
#' @return A numeric vector, matrix or data frame of the same dimension as
#' \code{data} is returned.
#' @note Please refer to the documentation for \code{box} for details about the
#' Box-Cox transformation in use.
#' @seealso \code{\link{box}}
#' @keywords math
#' @export 
rbox <- function(data, lambda) {
    if (length(lambda)>1 || lambda!=0) data <- sign(lambda*data+1)*(sign(lambda*data+1)*(lambda*data+1))^(1/lambda) else data <- exp(data)
    data
}


.ellipsePoints <- function(a,b, alpha = 0, loc = c(0,0), n = 100)
{
    ## Purpose: ellipse points,radially equispaced, given geometric par.s
    ## -------------------------------------------------------------------------
    ## Arguments: a, b : length of half axes in (x,y) direction
    ##            alpha: angle (in degrees) for rotation
    ##            loc  : center of ellipse
    ##            n    : number of points
    ## -------------------------------------------------------------------------
    ## Author: Martin Maechler, Date: 19 Mar 2002, 16:26
    ## modified by Kenneth to get rid of the precision problem met when there's a large difference in the length of the two axes

    small <- 0
    if (a/b  > 10 | b/a > 10) {
        ratio <- a/b
        b <- a
        if (round(alpha)==0) small <- 2 else small <- 1
    }

    B <- min(a,b)
    A <- max(a,b)
    ## B <= A
    d2 <- (A-B)*(A+B)                   #= A^2 - B^2
    
    phi <- 2*pi*seq(0,1, len = n)
    #phi<-t
    sp <- sin(phi)
    cp <- cos(phi)
    r <- a*b / sqrt(B^2 + d2 * sp^2)
    xy <- cbind(r,r) * cbind(cp, sp)
    ## xy are the ellipse points for alpha = 0 and loc = (0,0)
    al <- alpha * pi/180
    ca <- cos(al)
    sa <- sin(al)

 
    
    ######
    xy.new <- xy %*% rbind(c(ca, sa), c(-sa, ca))
    if (small==2) {
      #rotate after rescaling if appropriate
      #xy.new[,2]=xy.new[,2]/ratio
      xy.new <- xy.new %*% rbind(c(ca, sa), c(-sa/ratio, ca/ratio)) 
    }
    if (small==1) {
      #xy.new[,1]=xy.new[,1]/ratio
      xy.new <- xy.new %*% rbind(c(ca, sa), c(-sa/ratio, ca/ratio))
    }
    xy.new + cbind(rep(loc[1],n), rep(loc[2],n))
}

#' Scatterplot of Clustering Results
#' 
#' This method generates scatterplot revealing the cluster assignment, cluster
#' boundaries according to the specified percentile as well as supplemental
#' information like outliers or filtered observations.
#' 
#' 
#' @name plot,flowClust-method
#' @aliases plot
#' @docType methods
#' @param x Object returned from \code{\link{flowClust}}.
#' @param y missing
#' @param data A matrix, data frame of observations, or object of class
#' \code{flowFrame}. This is the object on which \code{flowClust} was
#' performed.
#' @param subset A numeric vector of length two indicating which two variables
#' are selected for the scatterplot.  Alternatively, a character vector
#' containing the names of the two variables is allowed if \code{x@varNames} is
#' not \code{NULL}.
#' @param ellipse A logical value indicating whether the cluster boundary is to
#' be drawn or not.  If \code{TRUE}, the boundary will be drawn according to
#' the level specified by \code{level} or \code{cutoff}.
#' @param show.outliers A logical value indicating whether outliers will be
#' explicitly shown or not.
#' @param show.rm A logical value indicating whether filtered observations will
#' be shown or not.
#' @param include A numeric vector specifying which clusters will be shown on
#' the plot.  By default, all clusters are included.
#' @param main Title of the plot.
#' @param grayscale A logical value specifying if a grayscale plot is desired.
#' This argument takes effect only if the default values of relevant graphical
#' arguments are taken.
#' @param col Color(s) of the plotting characters.  May specify a different
#' color for each cluster.
#' @param pch Plotting character(s) of the plotting characters.  May specify a
#' different character for each cluster.
#' @param cex Size of the plotting characters.  May specify a different size
#' for each cluster.
#' @param col.outliers Color of the plotting characters denoting outliers.
#' @param pch.outliers Plotting character(s) used to denote outliers.  May
#' specify a different character for each cluster.
#' @param cex.outliers Size of the plotting characters used to denote outliers.
#' May specify a different size for each cluster.
#' @param col.rm Color of the plotting characters denoting filtered
#' observations.
#' @param pch.rm Plotting character used to denote filtered observations.
#' @param cex.rm Size of the plotting character used to denote filtered
#' observations.
#' @param ecol Color(s) of the lines representing the cluster boundaries.  May
#' specify a different color for each cluster.
#' @param elty Line type(s) drawing the cluster boundaries.  May specify a
#' different line type for each cluster.
#' @param level,u.cutoff,z.cutoff These three optional arguments specify the
#' rule used to identify outliers.  By default, all of them are left
#' unspecified, meaning that the rule stated in \code{x@ruleOutliers} will be
#' taken.  Otherwise, these arguments will be passed to
#' \code{\link{ruleOutliers}}.
#' @param npoints The number of points used to draw each cluster boundary.
#' @param add A logical value.  If \code{TRUE}, add to the current plot.
#' @param \dots Further graphical parameters passed to the generic function
#' \code{plot}.
#' @note The cluster boundaries need not be elliptical since Box-Cox
#' transformation has been performed.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords graphs
#' @export 
#' @rdname plot.flowClust
setGeneric("plot", useAsDefault=plot)
#' @rdname plot.flowClust
#' @export 
setMethod("plot", signature(x="flowClust", y="missing"),
function(x, data, subset=c(1,2), ellipse=T, show.outliers=T, show.rm=F, include=1:(x@K), main=NULL, grayscale=F, col=(if (grayscale) gray(1/4) else 2:(length(include)+1)), pch=".", cex=0.6, col.outliers=gray(3/4), pch.outliers=".", cex.outliers=cex, col.rm=1, pch.rm=1, cex.rm=0.6, ecol=1, elty=1, level=NULL, u.cutoff=NULL, z.cutoff=NULL, npoints=100, add=F,...)
{
	if(ncol(x@mu)==1){
		hist(x,data,...)
		return(invisible(0));
	}
    if (is(data, "flowFrame")) data <- exprs(data)[,x@varNames]  else
    if (is(data, "matrix")) (if (length(x@varNames)>0) data <- as.matrix(data[,x@varNames]))  else
    if (is(data, "data.frame")) data <- as.matrix(data[,x@varNames])

    if (!is.numeric(subset)) subset <- match(subset, x@varNames)

    py <- ncol(data)
    data <- data[,subset]
		if(is.null(level)){
			level<-0.9
		}
    	label <- Map(x, rm.outliers=F)
    if (!add) plot(data, type="n", main=main, ...)  else title(main)
    flagFiltered <- is.na(label)

    # plot points with different colors/symbols corr. to cluster assignment
    col <- matrix(col, length(include))
    pch <- matrix(pch, length(include))
    cex <- matrix(cex, length(include))
    pch.outliers <- matrix(pch.outliers, length(include))
    cex.outliers <- matrix(cex.outliers, length(include))
    j <- 0
    if (!show.outliers) for (i in include)  points(data[!flagFiltered & label==i,], pch=pch[j <- j+1], col=col[j], cex=cex[j])  else {

        # plot outliers
        	if (!is.null(level) || !is.null(u.cutoff) || !is.null(z.cutoff)) ruleOutliers(x) <- list(level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff)
        	for (i in include) points(data[!flagFiltered & label==i & !x@flagOutliers,], pch=pch[j <- j+1], col=col[j], cex=cex[j])
        	j <- 0
        	for (i in include) points(data[!flagFiltered & label==i & x@flagOutliers,], pch=pch.outliers[j <- j+1], col=col.outliers, cex=cex.outliers[j])
    }

    # plot filtered points (from above or below)
    if (show.rm) points(data[flagFiltered,], pch=pch.rm, col=col.rm, cex=cex.rm)

    # plot ellipses
    if (ellipse) {
        ecol <- matrix(ecol, length(include))
        elty <- matrix(elty, length(include))

        if (all(x@nu!=Inf)) {
            if (x@ruleOutliers[1]==0) {     # 0 means quantile
				if(all(is.na(x@prior))){
                	cc <- py * qf(x@ruleOutliers[2], py, x@nu)
				}else{
					cc <- py * qf(x@ruleOutliers[2], py, x@nu)
				}
            }  else {     # 1 means u.cutoff
				if(all(is.na(x@prior))){
                	cc <- ((x@nu+py)/x@ruleOutliers[2] - x@nu)    
				}else{
					cc <- ((x@nu+py)/x@ruleOutliers[2] -x@nu)    
				}
            }
        }  else cc <- qchisq(x@ruleOutliers[2], py)

        j <- 0
        if (length(x@lambda)>0){
			if (any(x@lambda!=1)){
				lambda<-rep(x@lambda, length.out=x@K)
			}else{
			#WTF.. why set lambda to 0?
				lambda<-numeric(0)
			}
		}else{
			lambda<-numeric(0);
		}
        cc <- rep(cc, length.out=x@K)
        for (i in include) {
            eigenPair <- eigen(x@sigma[i,subset,subset])
            l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
            l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
            angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) * 180/pi

           if (length(lambda)>0&any(lambda!=1)) {
                points(rbox(.ellipsePoints(a=l1[i], b=l2[i], alpha=angle, loc=x@mu[i,subset], n=npoints), lambda[i]), type="l", lty=elty[j <- j+1], col=ecol[j])
            } else {
                points(.ellipsePoints(a=l1[i], b=l2[i], alpha=angle, loc=x@mu[i,subset], n=npoints), type="l", lty=elty[j <- j+1], col=ecol[j])
            }
        }  
    }

}
)

#' @rdname plot.flowClust
setMethod("plot", signature(x="flowClustList", y="missing"),
function(x, data, subset=c(1,2), ellipse=T, show.outliers=T, show.rm=F, include=1:(x@K), main=NULL, grayscale=F, col=(if (grayscale) gray(1/4) else 2:(length(include)+1)), pch=".", cex=0.6, col.outliers=gray(3/4), pch.outliers=".", cex.outliers=cex, col.rm=1, pch.rm=1, cex.rm=0.6, ecol=1, elty=1, level=NULL, u.cutoff=NULL, z.cutoff=NULL, npoints=501, add=F, ...)
{
    x <- as(x, "flowClust")
    selectMethod("plot", signature(x="flowClust", y="missing"))(x=x, data=data, subset=subset, ellipse=ellipse, show.outliers=show.outliers, show.rm=show.rm, include=include, main=main, grayscale=grayscale, col=col, pch=pch, cex=cex, col.outliers=col.outliers, pch.outliers=pch.outliers, cex.outliers=cex.outliers, col.rm=col.rm, pch.rm=pch.rm, cex.rm=cex.rm, ecol=ecol, elty=elty, level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff, npoints=npoints, add=add, ...)
}
)



# to compute the density of a multivariate t distribution with Box-Cox transformation


#' Density of the Multivariate t Distribution with Box-Cox Tranformation
#' 
#' This function computes the densities at the inputted points of the
#' multivariate \eqn{t} distribution with Box-Cox transformation.
#' 
#' 
#' @param x A matrix or data frame of size \eqn{N \times P}{N x P}, where
#' \eqn{N} is the number of observations and \eqn{P} is the dimension.  Each
#' row corresponds to one observation.
#' @param mu A numeric vector of length \eqn{P} specifying the mean.
#' @param sigma A matrix of size \eqn{P \times P}{P x P} specifying the
#' covariance matrix.
#' @param nu The degrees of freedom used for the \eqn{t} distribution.  If
#' \code{nu=Inf}, Gaussian distribution will be used.
#' @param lambda The Box-Cox transformation parameter.  If missing, the
#' conventional \eqn{t} distribution without transformation will be used.
#' @param log A logical value.  If \code{TRUE} then the logarithm of the
#' densities is returned.
#' @return A list with the following components: \item{value}{A vector of
#' length \eqn{N} containing the density values.} \item{md}{A vector of length
#' \eqn{N} containing the Mahalanobis distances.}
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @keywords distribution
#' @export 
dmvt <- function(x, mu, sigma, nu, lambda, log=FALSE) 
{
    if (is.vector(x) && length(x)==length(mu)) x <- matrix(x,1) else x <- as.matrix(x)
    p <- ncol(x)

   if (!missing(lambda)){
	if(lambda!=1){ 
		tx <- box(x, lambda)
	} else 
	{ 
	    tx <- x
	}
    } else 
	{
		tx <- x
   	} 	

    M <- mahalanobis(tx, mu, sigma)
    if (nu != Inf) value <- lgamma((nu+p)/2) - 1/2 * determinant(as.matrix(sigma), logarithm=T)$modulus[1] - p/2 * log(pi*nu) - lgamma(nu/2) - (nu+p)/2 * log(1+M/nu) else value <- -p/2 * log(2*pi) - 1/2 * determinant(as.matrix(sigma), logarithm=T)$modulus[1] - 1/2 * M

    # Jacobian of Box-Cox transformation
    # We ignore the Jacobian if no value for 'lambda' is specified.
    # We also ignore the Jacobian if 'lambda = 1', which corresponds to no
    # transformation. We do this to bypass the 'log(abs(x)))', which yields
    # '-Inf' whenever 'x = 0'. This, of course, is unintended when no
    # transformation is requested.
    if (!missing(lambda) && lambda != 1) {
      value <- value + (lambda-1) * rowSums(log(abs(x)))
    }
    if (log==F) value <- exp(value)
    list(value=value, md=M)
}


# to compute the density of a multivariate t mixture distribution with Box-Cox transformation


#' Density of the Multivariate t Mixture Distribution with Box-Cox
#' Tranformation
#' 
#' This function computes the densities at the inputted points of the
#' multivariate \eqn{t} mixture distribution with Box-Cox transformation.
#' 
#' 
#' @param x A matrix or data frame of size \eqn{N \times P}{N x P}, where
#' \eqn{N} is the number of observations and \eqn{P} is the dimension.  Each
#' row corresponds to one observation.
#' @param w A numeric vector of length \eqn{K} containing the cluster
#' proportions.
#' @param mu A matrix of size \eqn{K \times P}{K x P} containing the \eqn{K}
#' mean vectors.
#' @param sigma An array of size \eqn{K \times P \times P}{K x P x P}
#' containing the \eqn{K} covariance matrices.
#' @param nu A numeric vector of length \eqn{K} containing the degrees of
#' freedom used for the \eqn{t} distribution. If only one value is specified
#' for \code{nu}, then it is used for all \eqn{K} clusters. If \code{nu=Inf},
#' Gaussian distribution will be used.
#' @param lambda The Box-Cox transformation parameter.  If missing, the
#' conventional \eqn{t} distribution without transformation will be used.
#' @param object An optional argument.  If provided, it's an object returned
#' from \code{\link{flowClust}}, and the previous arguments will be assigned
#' values from the corresponding slots of \code{object}.
#' @param subset An optional argument.  If provided, it's a numeric vector
#' indicating which variables are selected for computing the densities.  If
#' \code{object} is provided and \code{object@varNames} is not \code{NULL},
#' then a character vector containing the names of the variables is allowed.
#' @param include An optional argument.  If provided, it's a numeric vector
#' specifying which clusters are included for computing the densities.
#' @param log A logical value.  If \code{TRUE} then the logarithm of the
#' densities is returned.
#' @return A vector of length \eqn{N} containing the density values.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @keywords distribution
#' @export 
dmvtmix <- function(x, w, mu, sigma, nu, lambda, object, subset, include, log=FALSE) 
{
    if (!missing(object)) {
        w <- object@w
        mu <- object@mu
        sigma <- object@sigma
        nu <- object@nu
        if (length(object@lambda)>0) 
            lambda <- object@lambda
        if (!missing(subset) && !is.numeric(subset)) 
            subset <- match(subset, object@varNames)
    }

    K <- length(w)
    if (K==1) {
        mu <- matrix(mu, 1)
        sigma <- array(sigma, c(1, ncol(mu), ncol(mu)))
    } else if (length(mu)==K) {
        mu <- matrix(mu, K, 1)
        sigma <- array(sigma, c(K, 1, 1))
    }

    # If only one value of 'nu' is specified, then it is repeated for each
    # mixture component. Otherwise, we allow the user to specify a different
    # value of 'nu' for each population. In the case that the number of values of
    # 'nu' is not the same as 'K', we throw an error.
    if (length(nu) == 1) {
      nu <- rep(nu, K)
    } else if (length(nu) != K) {
      stop("The number of values in 'nu' must be 1 or K.")
    }
    if (!missing(lambda)) 
        lambda <- rep(lambda, K)

    value <- 0
    if (missing(subset)) 
        subset <- 1:ncol(mu)
    if (missing(include)) 
        include <- 1:K
    sumw <- sum(w[include])
    for (k in include) {
        if (missing(lambda)) value <- value + w[k]/sumw * dmvt(x, mu[k,subset], sigma[k, subset, subset], nu[k])$value else value <- value + w[k]/sumw * dmvt(x, mu[k,subset], sigma[k, subset, subset], nu[k], lambda[k])$value
    }
    if (log) 
        value <- log(value)
    value
}



#if(!isGeneric("density")) setGeneric("density", useAsDefault=density)

#' Grid of Density Values for the Fitted t Mixture Model with Box-Cox
#' Transformation
#' 
#' This method constructs the \code{flowDens} object which is used to generate
#' a contour or image plot.
#' 
#' The \code{flowDens} object returned is to be passed to the \code{plot}
#' method for generating a contour or image plot.
#' 
#' @name density,flowClust-method
#' @aliases density-method density,flowClust-method
#' density,flowClustList-method flowDens-class density.flowClust
#' @docType methods
#' @param x Object returned from \code{\link{flowClust}} or from running
#' \code{\link[=filter.flowFrame]{filter}} on a \code{flowFrame} object.
#' @param data A matrix, data frame of observations, or object of class
#' \code{flowFrame}.  This is the object on which \code{flowClust} or
#' \code{filter} was performed.  If this argument is not specified, the grid
#' square upon which densities will be computed must be provided (through
#' arguments \code{from} and \code{to}).
#' @param subset A numeric vector of length two indicating which two variables
#' are selected for the scatterplot.  Alternatively, a character vector
#' containing the names of the two variables is allowed if \code{x@varNames} is
#' not \code{NULL}.
#' @param include A numeric vector specifying which clusters are included to
#' compute the density values.  By default, all clusters are included.
#' @param npoints A numeric vector of size two specifying the number of grid
#' points in \eqn{x} (horizontal) and \eqn{y} (vertical) directions
#' respectively.
#' @param from A numeric vector of size two specifying the coordinates of the
#' lower left point of the grid square.  Note that, if this (and \code{to}) is
#' not specified, \code{data} must be provided such that the range in the two
#' variables (dimensions) selected will be used to define the grid square.
#' @param to A numeric vector of size two specifying the co-ordinates of the
#' upper right point of the grid square.
#' @return An object of class \code{flowDens} containing the following slots is
#' constructed: \item{dx}{A numeric vector of length \code{npoints[1]}; the
#' \eqn{x}-coordinates of the grid points.} \item{dy}{A numeric vector of
#' length \code{npoints[2]}; the \eqn{y}-coordinates of the grid points.}
#' \item{value}{A matrix of size \code{npoints[1]} \eqn{\times}{x}
#' \code{npoints[2]}; the density values at the grid points.}
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link[=plot.flowDens]{plot}}, \code{\link{flowClust}}
#' @keywords graphs
#' @export 
#' @rdname density
#' @importFrom BiocGenerics density
#' @export 
setMethod("density", signature(x="flowClust"),
function(x, data=NULL, subset=c(1,2), include=1:(x@K), npoints=c(100,100), from=NULL, to=NULL)
{
    if (is(data, "flowFrame")) data <- exprs(data)[,x@varNames]  else
    if (is(data, "matrix")) (if (length(x@varNames)>0) data <- as.matrix(data[,x@varNames]))  else
    if (is(data, "data.frame")) data <- as.matrix(data[,x@varNames])

    if (length(colnames(data))==0) varNames <- NULL  else
    {
        if (!is.numeric(subset)) 
        {
            varNames <- subset
            subset <- match(subset, x@varNames)
        }
        else
        {
            varNames <- x@varNames[subset]
        }
    }

    trange <- (if (!is.null(data)) range(data[,subset[1]]) else c(from[1], to[1]))
    dx <- seq(from=min(trange), to=max(trange), by=abs(diff(trange))/(npoints[1]-1))
#    dx <- grid1(npoints[1], range=(if (!is.null(data)) range(data[,subset[1]]) else c(from[1], to[1])))
    trange <- (if (!is.null(data)) range(data[,subset[2]]) else c(from[2], to[2]))
    dy <- seq(from=min(trange), to=max(trange), by=abs(diff(trange))/(npoints[2]-1))
#    dy <- grid1(npoints[2], range=(if (!is.null(data)) range(data[,subset[2]]) else c(from[2], to[2])))
    xy <- cbind(rep(dx,length(dy)), rep(dy,each=length(dx)))
#    xy <- grid2(dx,dy)

    value <- 0
    nu <- rep(x@nu, length.out=x@K)
    if (length(x@lambda)>0&x@lambda!=1) {
        lambda <- rep(x@lambda, x@K)
        for (k in include) {
            xyTrans <- (apply(xy,2,sign)*apply(xy,2,abs)^lambda[k] - 1) / lambda[k]
            value <- value + x@w[k] * dmvt(xyTrans, x@mu[k,subset], x@sigma[k,subset,subset], nu[k], log=F)$value * abs(xy[,1])^(lambda[k]-1) * abs(xy[,2])^(lambda[k]-1)
        }
        value <- value / sum(x@w[include])
    } else {
        for (k in include) value <- value + x@w[k] * dmvt(xy, x@mu[k,subset], x@sigma[k,subset,subset], nu[k], log=F)$value
        value <- value / sum(x@w[include])
    }
    value <- matrix(value, length(dx), length(dy))

    dx <- matrix(dx)
    colnames(dx) <- varNames[1]
    dy <- matrix(dy)
    colnames(dy) <- varNames[2] 
    new("flowDens",dx=dx,dy=dy,value=value)
}
)

#' @rdname density
setMethod("density", signature(x="flowClustList"),
function(x, data=NULL, subset=c(1,2), include=1:(x@K), npoints=c(100,100), from=NULL, to=NULL)
{
    x <- as(x, "flowClust")
    callGeneric()
}
)


#' Contour or Image Plot of Clustering Results
#' 
#' This method makes use of the \code{flowDens} object returned by
#' \code{\link[=density.flowClust]{density}} to generate a contour or image
#' plot.
#' 
#' 
#' @name plot,flowDens-method
#' @aliases plot,flowDens,missing-method plot,flowDens-method plot.flowDens
#' @docType methods
#' @param x The \code{flowDens} object returned from
#' \code{\link[=density.flowClust]{density}}.
#' @param type Either \code{"contour"} or \code{"image"} to specify the type of
#' plot desired.
#' @param nlevels An integer to specify the number of contour levels or colors
#' shown in the plot.
#' @param scale If \code{"log"}, the logarithm of the density values will be
#' used to generate the plot; similar interpretation holds for \code{"sqrt"}.
#' The use of a \code{log} or \code{sqrt} elicits more information about low
#' density regions.
#' @param color A string containing the name of the function used to generate
#' the desired list of colors.
#' @param xlab,ylab Labels for the \eqn{x}- and \eqn{y}-axes respectively.
#' @param \dots Other arguments to be passed to \code{contour} or \code{image},
#' for example, \code{drawlabels} and \code{add}.  Once an image plot is
#' generated, users may impose a contour plot on it by calling this function
#' with an additional argument \code{add=TRUE}.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}, \code{\link[=density.flowClust]{density}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords graphs
#' @rdname plot.flowDens
setMethod("plot", signature(x="flowDens", y="missing"),
function(x, type=c("contour", "image"), nlevels=30, scale=c("raw", "log", "sqrt"), color=c("rainbow", "heat.colors", "terrain.colors", "topo.colors", "cm.colors", "gray"), xlab=colnames(x@dx), ylab=colnames(x@dy), ...)
{
    if (scale[1]=="log") z <- log(x@value)  else if (scale[1]=="sqrt") z <- sqrt(x@value)  else z <- x@value
    if (type[1]=="contour") contour(x=x@dx, y=x@dy, z=z, nlevels=nlevels, xlab=xlab, ylab=ylab, ...)  else {
        if (color[1]=="heat.colors") color <- heat.colors(nlevels) else
        if (color[1]=="rainbow") color <- rainbow(nlevels) else
        if (color[1]=="terrain.colors") color <- terrain.colors(nlevels) else
        if (color[1]=="topo.colors") color <- topo.colors(nlevels) else
        if (color[1]=="cm.colors") color <- cm.colors(nlevels) else
        if (color[1]=="gray") color <- gray(seq(0,1,length.out=nlevels))
        image(x=x@dx, y=x@dy, z=z, col=color, xlab=xlab, ylab=ylab, ...)
    }
}
)

