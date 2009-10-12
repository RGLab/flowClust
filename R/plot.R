# to perform box-cox transformation (multivariate)
box <- function(data, lambda) {
    if (length(lambda)>1 || lambda!=0) data <- (sign(data)*abs(data)^lambda-1)/lambda else data <- log(data)
    data
}


# to perform reverse box-cox transformation (multivariate)
rbox <- function(data, lambda) {
    if (length(lambda)>1 || lambda!=0) data <- sign(lambda*data+1)*(sign(lambda*data+1)*(lambda*data+1))^(1/lambda) else data <- exp(data)
    data
}


.ellipsePoints <- function(a,b, alpha = 0, loc = c(0,0), n = 501)
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
    if (a/b > 1000) {
        ratio <- a/b
        b <- a
        if (round(alpha)==0) small <- 2 else small <- 1
    }

    B <- min(a,b)
    A <- max(a,b)
    ## B <= A
    d2 <- (A-B)*(A+B)                   #= A^2 - B^2
    phi <- 2*pi*seq(0,1, len = n)
    sp <- sin(phi)
    cp <- cos(phi)
    r <- a*b / sqrt(B^2 + d2 * sp^2)
    xy <- r * cbind(cp, sp)
    ## xy are the ellipse points for alpha = 0 and loc = (0,0)
    al <- alpha * pi/180
    ca <- cos(al)
    sa <- sin(al)

    xy.new <- xy %*% rbind(c(ca, sa), c(-sa, ca))
    if (small==2) xy.new[,2]=xy.new[,2]/ratio
    if (small==1) xy.new[,1]=xy.new[,1]/ratio
    xy.new + cbind(rep(loc[1],n), rep(loc[2],n))
}



setMethod("plot", signature(x="flowClust", y="missing"),
function(x, data, subset=c(1,2), ellipse=T, show.outliers=T, show.rm=F, include=1:(x@K), main=NULL, grayscale=F, col=(if (grayscale) gray(1/4) else 2:(length(include)+1)), pch=".", cex=0.6, col.outliers=gray(3/4), pch.outliers=".", cex.outliers=cex, col.rm=1, pch.rm=1, cex.rm=0.6, ecol=1, elty=1, level=NULL, u.cutoff=NULL, z.cutoff=NULL, npoints=501, add=F, ...)
{
    if (is(data, "flowFrame")) data <- exprs(data)[,x@varNames]  else
    if (is(data, "matrix")) (if (length(x@varNames)>0) data <- as.matrix(data[,x@varNames]))  else
    if (is(data, "data.frame")) data <- as.matrix(data[,x@varNames])

    if (!is.numeric(subset)) subset <- match(subset, x@varNames)

    py <- ncol(data)
    data <- data[,subset]
    label <- map(x@z)
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
                cc <- py * qf(x@ruleOutliers[2], py, x@nu)
            }  else {     # 1 means u.cutoff
                cc <- ((x@nu+py)/x@ruleOutliers[2] - x@nu)    
            }
        }  else cc <- qchisq(x@ruleOutliers[2], py)

        j <- 0
        lambda <- if (length(x@lambda)>0) rep(x@lambda, length.out=x@K) else numeric(0)
        cc <- rep(cc, length.out=x@K)
        for (i in include) {
            eigenPair <- eigen(x@sigma[i,subset,subset])
            l1 <- sqrt(eigenPair$values[1]) * sqrt(cc)
            l2 <- sqrt(eigenPair$values[2]) * sqrt(cc)
            angle <- atan(eigenPair$vectors[2,1] / eigenPair$vectors[1,1]) * 180/pi

            if (length(lambda)>0) {
                points(rbox(.ellipsePoints(a=l1[i], b=l2[i], alpha=angle, loc=x@mu[i,subset], n=npoints), lambda[i]), type="l", lty=elty[j <- j+1], col=ecol[j])
            } else {
                points(.ellipsePoints(a=l1[i], b=l2[i], alpha=angle, loc=x@mu[i,subset], n=npoints), type="l", lty=elty[j <- j+1], col=ecol[j])
            }
        }  
    }

}
)


setMethod("plot", signature(x="flowClustList", y="missing"),
function(x, data, subset=c(1,2), ellipse=T, show.outliers=T, show.rm=F, include=1:(x@K), main=NULL, grayscale=F, col=(if (grayscale) gray(1/4) else 2:(length(include)+1)), pch=".", cex=0.6, col.outliers=gray(3/4), pch.outliers=".", cex.outliers=cex, col.rm=1, pch.rm=1, cex.rm=0.6, ecol=1, elty=1, level=NULL, u.cutoff=NULL, z.cutoff=NULL, npoints=501, add=F, ...)
{
    x <- as(x, "flowClust")
    selectMethod("plot", signature(x="flowClust", y="missing"))(x=x, data=data, subset=subset, ellipse=ellipse, show.outliers=show.outliers, show.rm=show.rm, include=include, main=main, grayscale=grayscale, col=col, pch=pch, cex=cex, col.outliers=col.outliers, pch.outliers=pch.outliers, cex.outliers=cex.outliers, col.rm=col.rm, pch.rm=pch.rm, cex.rm=cex.rm, ecol=ecol, elty=elty, level=level, u.cutoff=u.cutoff, z.cutoff=z.cutoff, npoints=npoints, add=add, ...)
}
)



# to compute the density of a multivariate t distribution with Box-Cox transformation
dmvt <- function(x, mu, sigma, nu, lambda, log=FALSE) 
{
    if (is.vector(x) && length(x)==length(mu)) x <- matrix(x,1) else x <- as.matrix(x)
    p <- ncol(x)

    if (!missing(lambda)) tx <- box(x, lambda) else tx <- x

    M <- mahalanobis(tx, mu, sigma)
    if (nu != Inf) value <- lgamma((nu+p)/2) - 1/2 * determinant(as.matrix(sigma), log=T)$modulus[1] - p/2 * log(pi*nu) - lgamma(nu/2) - (nu+p)/2 * log(1+M/nu) else value <- -p/2 * log(2*pi) - 1/2 * determinant(as.matrix(sigma), log=T)$modulus[1] - 1/2 * M
    if (!missing(lambda)) value <- value + (lambda-1) * rowSums(log(abs(x)))
    if (log==F) value <- exp(value)
    list(value=value, md=M)
}


# to compute the density of a multivariate t mixture distribution with Box-Cox transformation
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

    nu <- rep(nu, K)
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



if(!isGeneric("density")) setGeneric("density", useAsDefault=density)


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

    dx <- grid1(npoints[1], range=(if (!is.null(data)) range(data[,subset[1]]) else c(from[1], to[1])))
    dy <- grid1(npoints[2], range=(if (!is.null(data)) range(data[,subset[2]]) else c(from[2], to[2])))
    xy <- grid2(dx,dy)

    value <- 0
    nu <- rep(x@nu, length.out=x@K)
    if (length(x@lambda)>0) {
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


setMethod("density", signature(x="flowClustList"),
function(x, data=NULL, subset=c(1,2), include=1:(x@K), npoints=c(100,100), from=NULL, to=NULL)
{
    x <- as(x, "flowClust")
    callGeneric()
}
)



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



if(!isGeneric("hist")) setGeneric("hist",useAsDefault=hist)


setMethod("hist", signature(x="flowClust"),
function(x, data=NULL, subset=1, include=1:(x@K), histogram=TRUE, labels=TRUE, xlim=NULL, ylim=NULL, xlab=(if (is.numeric(subset)) NULL else subset), ylab="Density", main=NULL, breaks=50, col=NULL, pch=20, cex=0.6, ...)
{
    den <- function(y) {
        value <- 0
        nu <- rep(x@nu, length.out=x@K)
        if (length(x@lambda)>0) {
            lambda <- rep(x@lambda, length.out=x@K)
            for (k in include) {
                yTrans <- (sign(y)*abs(y)^lambda[k] - 1) / lambda[k]
                value <- value + x@w[k] * dmvt(yTrans, x@mu[k,subset], x@sigma[k,subset,subset], nu[k], log=F)$value * abs(y)^(lambda[k]-1)
            }
        } else {
            for (k in include) value <- value + x@w[k] * dmvt(y, x@mu[k,subset], x@sigma[k,subset,subset], nu[k], log=F)$value
        }
        value <- value / sum(x@w[include])
        value
    }

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
        ymax <- max(den(tseq[tseq!=0]))
    }

    # look for highest point in histogram
    if (histogram) {
        data2 <- data[!is.na(x@flagOutliers) & is.element(map(x@z), include)]
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
    
    if (histogram) hist(data2, breaks=tbreaks, freq=F, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...)

    curve(den, add=histogram, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, main=main, ...)

    if (labels) {
        if (is.null(col)) {
            if (length(include)<=4) col <- c("red", "blue", "green", "black")  else col <- 2:(length(include)+1)
        } else col<-matrix(col, length(include))
        j <- 0
        for (k in include) stripchart(data[map(x@z)==k], add=T, at=ymin - (ylim[2]-ymin)/100*(j<-j+1), pch=pch, cex=cex, col=col[j])
    }
}
)


setMethod("hist", signature(x="flowClustList"),
function(x, data=NULL, subset=1, include=1:(x@K), histogram=T, labels=T, xlim=NULL, ylim=NULL, xlab=(if (is.numeric(subset)) NULL else subset), ylab="Density", main=NULL, breaks=50, col=NULL, pch=20, cex=0.6, ...)
{
    x <- as(x, "flowClust")
    callGeneric()
}
)
