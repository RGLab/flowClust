#' Robust Model-based Clustering for Flow Cytometry
#' 
#' This function performs automated clustering for identifying cell populations
#' in flow cytometry data.  The approach is based on the tmixture model
#' with the Box-Cox transformation, which provides a unified framework to
#' handle outlier identification and data transformation simultaneously.
#' 
#' Estimation of the unknown parameters (including the Box-Cox parameter) is
#' done via an Expectation-Maximization (EM) algorithm.  At each EM iteration,
#' Brent's algorithm is used to find the optimal value of the Box-Cox
#' transformation parameter.  Conditional on the transformation parameter, all
#' other estimates can be obtained in closed form.  Please refer to Lo et al.
#' (2008) for more details.
#' 
#' The \pkg{flowClust} package makes extensive use of the GSL as well as BLAS.
#' If an optimized BLAS library is provided when compiling the package, the
#' \pkg{flowClust} package will be able to run multi-threaded processes.
#' 
#' Various operations have been defined for the object returned from
#' \code{\link{flowClust}}.  
#' 
#' In addition, to facilitate the integration with the \pkg{flowCore} package
#' for processing flow cytometry data, the \code{flowClust} operation can be
#' done through a method pair (\code{\link{tmixFilter}} and
#' \code{\link[=tmixFilter]{filter}}) such that various methods defined in
#' \pkg{flowCore} can be applied on the object created from the filtering
#' operation.
#' 
#' @name fowClust
#' 
#' @param x A numeric vector, matrix, data frame of observations, or object of
#' class \code{flowFrame}.  Rows correspond to observations and columns
#' correspond to variables.
#' @param expName A character string giving the name of the experiment.
#' @param varNames A character vector specifying the variables (columns) to be
#' included in clustering.  When it is left unspecified, all the variables will
#' be used.
#' @param K An integer vector indicating the numbers of clusters.
#' @param nu The degrees of freedom used for the \eqn{t} distribution.  Default
#' is 4.  If \code{nu=Inf}, Gaussian distribution will be used.
#' @param lambda The initial transformation to be applied to the data.
#' @param trans A numeric indicating whether the Box-Cox transformation
#' parameter is estimated from the data.  May take 0 (no estimation), 1
#' (estimation, default) or 2 (cluster-specific estimation).
#' @param min.count An integer specifying the threshold count for filtering
#' data points from below.  The default is 10, meaning that if 10 or more data
#' points are smaller than or equal to \code{min}, they will be excluded from
#' the analysis.  If \code{min} is \code{NULL}, then the minimum of data as per
#' each variable will be used.  To suppress filtering, set it as -1.
#' @param max.count An integer specifying the threshold count for filtering
#' data points from above.  Interpretation is similar to that of
#' \code{min.count}.
#' @param min The lower boundary set for data filtering.  Note that it is a
#' vector of length equal to the number of variables (columns), implying that a
#' different value can be set as per each variable.
#' @param max The upper boundary set for data filtering.  Interpretation is
#' similar to that of \code{min}.
#' @param randomStart A numeric value indicating how many times a random
#' parition of the data is generated for initialization.  The default is 0,
#' meaning that a deterministic partition based on kmeans clustering is used. A
#' value of 10 means random partitions of the data will be generated, each of
#' which is followed by a short EM run.  The partition leading to the highest
#' likelihood value will be adopted to be the initial partition for the
#' eventual long EM run.
#' @param prior The specification of the prior. Used if usePrior="yes"
#' @param usePrior Argument specifying whether or not the prior will be used.
#' Can be "yes","no","vague". A vague prior will be automatically specified if
#' usePrior="vague"
#' @param criterion A character string stating the criterion used to choose the
#' best model.  May take either \code{"BIC"} or \code{"ICL"}.  This argument is
#' only relevant when \code{length(K)>1}. Default is "BIC".
#' @param ... other arguments: B: The maximum number of EM iterations.Default
#' is 500.
#' 
#' tol: The tolerance used to assess the convergence of the EM. default is
#' 1e-5.
#' 
#' nu.est: A numeric indicating whether \code{nu} is to be estimated or not.
#' May take 0 (no estimation, default), 1 (estimation) or 2 (cluster-specific
#' estimation). Default is 0.
#' 
#' level: A numeric value between 0 and 1 specifying the threshold quantile
#' level used to call a point an outlier.  The default is 0.9, meaning that any
#' point outside the 90\% quantile region will be called an outlier.
#' 
#' u.cutoff: Another criterion used to identify outliers.  If this is
#' \code{NULL}, which is default, then \code{level} will be used.  Otherwise,
#' this specifies the threshold (e.g., 0.5) for \eqn{u}, a quantity used to
#' measure the degree of \dQuote{outlyingness} based on the Mahalanobis
#' distance.  Please refer to Lo et al. (2008) for more details.
#' 
#' z.cutoff: A numeric value between 0 and 1 underlying a criterion which may
#' be used together with \code{level}/\code{u.cutoff} to identify outliers.  A
#' point with the probability of assignment \eqn{z} (i.e., the posterior
#' probability that a data point belongs to the cluster assigned) smaller than
#' \code{z.cutoff} will be called an outlier.  The default is 0, meaning that
#' assignment will be made no matter how small the associated probability is,
#' and outliers will be identified solely based on the rule set by \code{level}
#' or \code{cutoff}.
#' 
#' B.init: The maximum number of EM iterations following each random partition
#' in random initialization. Default is the same as B.
#' 
#' tol.init: The tolerance used as the stopping criterion for the short EM runs
#' in random initialization. Default is 1e-2.
#' 
#' seed: An integer giving the seed number used when
#' \code{randomStart>0}.Default is 1.
#' 
#' control: An argument reserved for internal use.
#' @return If \code{K} is of length 1, the function returns an object of class
#' \code{flowClust} containing the following slots, where \eqn{K} is the number
#' of clusters, \eqn{N} is the number of observations and \eqn{P} is the number
#' of variables: \item{expName}{Content of the \code{expName} argument.}
#' \item{varNames}{Content of the \code{varNames} argument if provided;
#' generated if available otherwise.} \item{K}{An integer showing the number of
#' clusters.} \item{w}{A vector of length \eqn{K}, containing the estimates of
#' the \eqn{K} cluster proportions.} \item{mu}{A matrix of size \eqn{K \times
#' P}{K x P}, containing the estimates of the \eqn{K} mean vectors.}
#' \item{sigma}{An array of dimension \eqn{K \times P \times P}{K x P x P},
#' containing the estimates of the \eqn{K} covariance matrices.}
#' \item{lambda}{The Box-Cox transformation parameter estimate.} \item{nu}{The
#' degrees of freedom for the \eqn{t} distribution.} \item{z}{A matrix of size
#' \eqn{N \times K}{N x K}, containing the posterior probabilities of cluster
#' memberships.  The probabilities in each row sum up to one.} \item{u}{A
#' matrix of size \eqn{N \times K}{N x K}, containing the \dQuote{weights} (the
#' contribution for computing cluster mean and covariance matrix) of each data
#' point in each cluster.  Since this quantity decreases monotonically with the
#' Mahalanobis distance, it can also be interpreted as the level of
#' \dQuote{outlyingness} of a data point.  Note that, when \code{nu=Inf}, this
#' slot is used to store the Mahalanobis distances instead.} \item{label}{A
#' vector of size \eqn{N}, showing the cluster membership according to the
#' initial partition (i.e., hierarchical clustering if \code{randomStart=0} or
#' random partitioning if \code{randomStart>0}).  Filtered observations will be
#' labelled as \code{NA}.  Unassigned observations (which may occur since only
#' 1500 observations at maximum are taken for hierarchical clustering) will be
#' labelled as 0.} \item{uncertainty}{A vector of size \eqn{N}, containing the
#' uncertainty about the cluster assignment.  Uncertainty is defined as 1 minus
#' the posterior probability that a data point belongs to the cluster to which
#' it is assigned.} \item{ruleOutliers}{A numeric vector of size 3, storing the
#' rule used to call outliers.  The first element is 0 if the criterion is set
#' by the \code{level} argument, or 1 if it is set by \code{u.cutoff}.  The
#' second element copies the content of either the \code{level} or
#' \code{u.cutoff} argument.  The third element copies the content of the
#' \code{z.cutoff} argument.  For instance, if points are called outliers when
#' they lie outside the 90\% quantile region or have assignment probabilities
#' less than 0.5, then \code{ruleOutliers} is \code{c(0, 0.9, 0.5)}.  If points
#' are called outliers only if their \dQuote{weights} in the assigned clusters
#' are less than 0.5 regardless of the assignment probabilities, then
#' \code{ruleOutliers} becomes \code{c(1, 0.5, 0)}.} \item{flagOutliers}{A
#' logical vector of size \eqn{N}, showing whether each data point is called an
#' outlier or not based on the rule defined by \code{level}/\code{u.cutoff} and
#' \code{z.cutoff}.} \item{rm.min}{Number of points filtered from below.}
#' \item{rm.max}{Number of points filtered from above.} \item{logLike}{The
#' log-likelihood of the fitted mixture model.} \item{BIC}{The Bayesian
#' Information Criterion for the fitted mixture model.} \item{ICL}{The
#' Integrated Completed Likelihood for the fitted mixture model.} If \code{K}
#' has a length >1, the function returns an object of class
#' \code{flowClustList}.  Its data part is a list with the same length as
#' \code{K}, each element of which is a \code{flowClust} object corresponding
#' to a specific number of clusters.  In addition, the resultant
#' \code{flowClustList} object contains the following slots:\cr
#' 
#' \code{index} An integer giving the index of the list element corresponding
#' to the best model as selected by \code{criterion}.\cr \code{criterion} The
#' criterion used to choose the best model -- either \code{"BIC"} or
#' \code{"ICL"}.\cr
#' 
#' Note that when a \code{flowClustList} object is used in place of a
#' \code{flowClust} object, in most cases the list element corresponding to the
#' best model will be extracted and passed to the method/function call.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link[=summary.flowClust]{summary}},
#' \code{\link[=plot.flowClust]{plot}},
#' \code{\link[=density.flowClust]{density}},
#' \code{\link[=hist.flowClust]{hist}}, \code{\link{Subset}},
#' \code{\link{split}}, \code{\link{ruleOutliers}}, \code{\link{Map}},
#' \code{\link{SimulateMixture}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords cluster models
#' @examples
#' 
#' data(rituximab)
#' 
#' ### cluster the data using FSC.H and SSC.H
#' res1 <- flowClust(rituximab, varNames=c("FSC.H", "SSC.H"), K=1)
#' 
#' ### remove outliers before proceeding to the second stage
#' # %in% operator returns a logical vector indicating whether each
#' # of the observations lies within the cluster boundary or not
#' rituximab2 <- rituximab[rituximab %in% res1,]
#' # a shorthand for the above line
#' rituximab2 <- rituximab[res1,]
#' # this can also be done using the Subset method
#' rituximab2 <- Subset(rituximab, res1)
#' 
#' ### cluster the data using FL1.H and FL3.H (with 3 clusters)
#' res2 <- flowClust(rituximab2, varNames=c("FL1.H", "FL3.H"), K=3)
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
#' # transformed back to the original scale)
#' getEstimates(res2)$locations
#' 
#' ### demonstrate the use of various plotting methods
#' # a scatterplot
#' plot(res2, data=rituximab2, level=0.8)
#' plot(res2, data=rituximab2, level=0.8, include=c(1,2), grayscale=TRUE,
#'     pch.outliers=2)
#' # a contour / image plot
#' res2.den <- density(res2, data=rituximab2)
#' plot(res2.den)
#' plot(res2.den, scale="sqrt", drawlabels=FALSE)
#' plot(res2.den, type="image", nlevels=100)
#' plot(density(res2, include=c(1,2), from=c(0,0), to=c(400,600)))
#' # a histogram (1-D density) plot
#' hist(res2, data=rituximab2, subset="FL1.H")
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
#' flowClust(rituximab2, varNames=c("FL1.H", "FL3.H"), K=3, B=100, 
#'     min=c(0,0), max=c(400,800), level=0.95, z.cutoff=0.5)
#' @rdname flowClust
#' @import methods graph flowCore
#' @importFrom parallel mclapply
#' @useDynLib flowClust
#' @export 
flowClust<-function(x, expName="Flow Experiment", varNames=NULL, K
                    , nu=4, lambda=1,trans=1, min.count=10, max.count=10, min=NULL, max=NULL
                    ,  randomStart=0, prior=NULL,usePrior="no", criterion="BIC", ...)
{
	if (is(x, "flowFrame")) {
		if (length(varNames)==0) {
			y <- exprs(x)
			varNames <- colnames(y)
		}
		else {
			y <- as.matrix(exprs(x)[, varNames,drop=FALSE])
		}
	}
	else if (is(x, "matrix")) {
		if (length(varNames)==0) {
			y <- x
			if (length(colnames(x))==0)
				varNames <- "Not Available"
			else varNames <- colnames(x)
		}
		else {
			y <- as.matrix(x[, varNames,drop=FALSE])
		}
	}
	else if (is(x, "data.frame")) {
		if (length(varNames)==0) {
			y <- as.matrix(x)
			varNames <- colnames(x)
		}
		else {
			y <- as.matrix(x[, varNames,drop=FALSE])
		}
	}
	else if (is(x, "vector")) {
		y <- matrix(x)
		if (length(varNames)==0)
			varNames <- "Not Available"
	}
	else {
		stop(paste("Object ", as.character(x), " is not of class flowFrame / matrix / data frame!"))
	}

	# finding filtered observations
	rm.max <- rm.min <- rep(FALSE, nrow(y))
	if (max.count > -1) {
		if (is.null(max)[1])
			max <- apply(y, 2, max)
		for (k in 1:ncol(y))  if (sum(y[,k]>=max[k]) >= max.count)
				rm.max <- rm.max | (y[,k] >= max[k])
	}
	if (min.count > -1) {
		if (is.null(min)[1])
			min <- apply(y, 2, min)
		for (k in 1:ncol(y))  if (sum(y[,k]<=min[k]) >= min.count)
				rm.min <- rm.min | (y[,k] <= min[k])
	}
	if(min.count <= -1 && max.count <= -1)
	  include <- NULL
	else
	  include <- !rm.max & !rm.min

	usePrior=match.arg(as.character(usePrior),c("yes","no"))
	if(usePrior=="yes"){
		if(is.null(prior)){
			stop("You must specify a prior with usePrior=\"yes\"");
		}
		if (length(K)>1)
		{
			stop("You can only use a prior if the number of cluster is fixed!")
		}
		if(randomStart){
			message("randomStart>0 has no effect when using a prior. Labels are initialized from the prior.");
		}
		if(trans==2)
		{
			stop("You are using a prior with cluster specific transformations.\n  This is not recommended.")
		}
                if(any(is.infinite(c(nu,prior$nu0)))){
                    stop("If usePrior='yes', nu or nu0 may not be Inf");
                }
                if(lambda!=1&trans==1){
                    warning("Use of a prior with transformation estimation and lambda!=1 requires the prior means to be on the transformed scale.")
                }
		# TODO Add tests for validity of w0. Set a default. Same for oorder.
		if(!is.null(prior)&length(K)==1){
			#Check that the prior dimensions match the model being fit.
			mp<-ncol(prior$Mu0);
			mk<-nrow(prior$Mu0);
			ld<-dim(prior$Lambda0)
			lomega<-dim(prior$Omega0);
			#message(lomega);
			lnu0<-length(prior$nu0);
			py<-ncol(y);
			warn.o<-options("warn")
			options(warn=-1);
			if(any(is.null(prior$w0)))stop("w0 prior should be defined")
			if(length(prior$w0)!=K|(length(ld)!=3)|((lnu0!=1)&(lnu0!=mk))|(mp!=py)|(mk!=K)|(ld[1]!=K)|(ld[2]!=py)|(ld[3]!=py)|((any(lomega!=c(K,py,py)))&(any(lomega!=c(py,py))))){
				stop("Prior dimensions must match the number of clusters and the number of dimensions being analyzed.")
			}
			if(lnu0==1){
				#Extend nu0 to be a vector of length k. We allow cluster specific weight for covariance matrices.
				prior$nu0<-rep(prior$nu0,mk);
			}
			if(length(prior$w0)==1){
					#Extend w0 to be a vector of length k.
					prior$w0<-rep(prior$w0,mk);
			}
			options(warn=warn.o$warn)
		}
	}
	if (usePrior=="no")
	{
		if(!is.null(prior)){
			message("The prior specification has no effect when usePrior=",usePrior);
		}
		prior<-list(NA);
	}
	mc.cores <- getOption("mc.cores", 2L)
	if(mc.cores < 2||length(K) == 1)
	{
		message("Using the serial version of flowClust")
		# C version
		result<-lapply(as.list(1:length(K)),.flowClustK, y, expName=expName, varNames=varNames, K=K, criterion=criterion
        , nu=nu, lambda=lambda, trans=trans, min.count=min.count, max.count=max.count, min=min, max=max
        , randomStart=randomStart, include=include, rm.max, rm.min, prior,usePrior
        , ...)
	}
	else{
      
		# Split into nClust segReadsList
      # We solely rely on getOption("mc.cores",2L) to determine parallel cores.
      # and don't want to pass mc.cores explicitly because on windows, mclapply does not take mc.cores>1 
		result<-mclapply(as.list(1:length(K)),.flowClustK, y, expName=expName, varNames=varNames, K=K, criterion=criterion
        , nu=nu, lambda=lambda, trans=trans, min.count=min.count, max.count=max.count, min=min, max=max
        , randomStart=randomStart, include=include, rm.max, rm.min, prior,usePrior, mc.preschedule=FALSE
        , ...)
	}
	  
	# Simply return a flowClust object
	if (length(K)==1)
	{
		result[[1]]
	}
	# Create a list flowClustList
	else
	{
		result <- new("flowClustList", result, criterion=criterion)
		result@index <- which.max(criterion(result, criterion))
		result
	}
}

.flowClustK<-function(i, y, expName="Flow Experiment", varNames=NULL, K
                        , nu, lambda, trans, min.count, max.count, min, max, randomStart, include, rm.max,rm.min, prior,usePrior, criterion # default values set in flowClust API
                        , nu.est=0, B=500, tol=1e-5, level=0.9, u.cutoff=NULL, z.cutoff=0, B.init=B, tol.init=1e-2, seed=1, control=NULL
                        , nstart = 100 # passed to kmeans, only has effect on more than 1d clustering
                    )
{
	oorder<-1:K[i]
	.model<-1; #Tells the C code whether to run ECM with non-conjugate priors, or classic flowClust.'
	match.arg(usePrior,c("yes","no"));
	switch(usePrior,
			yes=priorFlag<-1,
			no=priorFlag<-0)
	if(!is.null(include))
	  y <- as.matrix(y[include,,drop=FALSE])
	ly <- nrow(y)
	py <- ncol(y)
	if (min(y)<=0 && lambda<=0)
		stop("lambda must be positive when data contain zero / negative values!")
	else if(usePrior=="yes")
	{
    #If we have a prior lambda.. use it to override the specified lambda.
    if(exists("lambda",prior)){
      if(!is.null(prior$lambda)){
        if(prior$lambda!=0){
          lambda<-prior$lambda
        }
      }
    }
		Mu0<-prior$Mu0
		Lambda0<-prior$Lambda0
		if(length(dim(prior$Omega0))==2){
			Omega0<-array(prior$Omega0,c(py,py,K[i]))
		}else{
			Omega0<-aperm(prior$Omega0,c(2:3,1));
		}

		for(j in 1:K[i]){
			if(!all(as.vector(Omega0[,,j])==0)){
				#message(j)
				#message(Omega0[,,j])
				Omega0[,,j]<-try(solve((Omega0[,,j])))
			}else{
				#message(j);
				#message(Omega0[,,j])
				Omega0[,,j]<-Omega0[,,j];
			}
		}
		nu0<-prior$nu0
		w0<-prior$w0;
		kappa0<-0.0
		Lambda0<-aperm(Lambda0,c(2:3,1))
		initprec<-array(0,c(py,py,K[i]));
		#cat("initprec dim: ",dim(initprec),"\n")
		#cat("Lambda0 dim:",dim(Lambda0),"\n");
		if(!all(as.vector(Lambda0==0))){
			for(j in 1:dim(Lambda0)[3]){
				#cat("Lambda0 value:",(Lambda0[,,j]),"\n")
				#cat("initprec value:",initprec[,,j],"\n")
				initprec[,,j]<-try(chol(solve(Lambda0[,,j])),silent=TRUE);
				#cat("initprec class:",class(initprec[,,]))
				#cat("solve initprec value:",initprec[,,j],"\n")
				if(inherits(initprec[,,j],"try-error")){
					##Should do something here to catch the error.
				}
			}
		}else{
			initprec<-array(var(box(y,lambda))/(K[i]^(2/py)),c(py,py,K[i]));
		}
		priorFlag<-1
		.model<-2;
	}
	else
	{
		# Non informative prior
		Mu0<-matrix(rep(colMeans(y),K[i]),K[i],py,byrow=TRUE)
		nu0<- rep(-py-1,K[i]);
		#Don't need Omega0, we'll use the conjugate prior code.
		Omega0<-rep(0,py*py);
		kappa0<-0;
		initprec<-rep(0,K[i]*py*py);
		#Lambda0, the prior covariance
		Lambda0<-rep(0,K[i]*py*py)
		w0<-rep(0,K[i]);
		.model<-1;
	}

# to determine the rule of calling outliers
	if (nu != Inf) {
		if (is.null(u.cutoff)) {
			if (!nu.est) {
				cc <- py * qf(level, py, nu)
				u.cutoff <- (nu + py) / (nu + cc)
			}
			else {
				u.cutoff <- level     # pass level to .C call to compute u.cutoff
			}
			ruleOutliers <- c(0, level, z.cutoff)     # 0 means quantile
		}
		else {
			ruleOutliers <- c(1, u.cutoff, z.cutoff)     # 1 means cutoff
		}
	}
	else {
		if (level != 1)
			q.cutoff <- qchisq(level, py)
		else
			q.cutoff <- -1     # -1 means no outlier identification
		ruleOutliers <- c(0, level, z.cutoff)
	}


	if (is.null(control$B.lambda)) control$B.lambda <- B    # BSolve=100
	if (is.null(control$B.brent)) control$B.brent <- 10000    # iterSolveMax=50
	if (is.null(control$tol.brent)) control$tol.brent <- 1e-5    # DiffSolve=1e-3
	if (is.null(control$xLow)) control$xLow <- 0.1    # xLow=.1
	if (is.null(control$xUp)) control$xUp <- 10    # xUp=1
	if (is.null(control$nuLow)) control$nuLow <- 2    # nuLow=2
	if (is.null(control$nuUp)) control$nuUp <- 100    # nuUp=30
	if (is.null(control$seed)) control$seed <- TRUE

	ind <- 0
	if (K[i]==1)
	{
		label <- rep(1, ly)
	}
	else if (!randomStart) #kmeans initialization
	{
		if (py==1)
		{
			q <- quantile(y, seq(from=0, to=1, by=1/K[i]))
			label <- rep(0, ly)
			q[1] <- q[1]-1
			for (k in 1:K[i]) label[y>q[k] & y<=q[k+1]] <- k
		}else{
			#label<-try(kmeans(y,Mu0)$cluster,silent=TRUE)
			#if(inherits(label,"try-error"))
			label<-try(kmeans(scale(y),K[i],nstart=nstart,iter.max=100)$cluster,silent=TRUE)
		}
	}
	if(priorFlag==0)
	{ # Initialization based on short EMs with random partitions if randomStart=TRUE
		if (control$seed) set.seed(seed)
		if (randomStart==1)
		{
			label <- sample(1:K[i], ly, replace=T)
		}
		else if(randomStart==0)
		{
			M<-0
			if(inherits(label,"try-error")){
				label <- sample(1:K[i], ly, replace=T)
			}
			ind<-0;
		}
		else
		{
			maxLabel <- vector("list",randomStart)
			maxLogLike <- rep(NA,randomStart)
			for (j in 1:randomStart)
			{
				label <- sample(1:K[i], ly, replace=TRUE)
				if (nu != Inf)
				{
					#cat(initprec,"\n");
					#ordering of the priors.. used to reorder the population names;
					obj <- try(.C("flowClust", as.double(t(y)), as.integer(ly),
									as.integer(py), as.integer(K[i]),
									w=rep(0,K[i]), mu=rep(0,K[i]*py),
									precision=initprec,
									lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
									nu=as.double(rep(nu,K[i])),
									z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
									as.integer(label), uncertainty=double(ly),
									as.double(rep(u.cutoff,K[i])), as.double(z.cutoff),
									flagOutliers=integer(ly), as.integer(B.init),
									as.double(tol.init), as.integer(trans),
									as.integer(nu.est), logLike=as.double(0),
									as.integer(control$B.lambda), as.integer(control$B.brent),
									as.double(control$tol.brent), as.double(control$xLow),
									as.double(control$xUp), as.double(control$nuLow),
									as.double(control$nuUp),
									mu0=as.double(t(Mu0)),
									as.double(kappa0),
									nu0=as.double(nu0),
									lambda0=as.double(Lambda0),
									omega0=as.double(Omega0),
									w0=as.double(w0),
									as.integer(.model),
									oorder=as.integer(oorder),
									PACKAGE="flowClust"))
									if (class(obj)=="try-error")
									{
										message("flowClust failed")
									}
				}
				else
				{
					obj <- try(.C("flowClustGaussian", as.double(t(y)), as.integer(ly),
									as.integer(py), as.integer(K[i]),
									w=rep(0,K[i]), mu=rep(0,K[i]*py),
									precision=initprec,
									lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
									z=rep(0, ly*K[i]), u=rep(0,ly*K[i]),
									as.integer(label), uncertainty=double(ly),
									as.double(q.cutoff), as.double(z.cutoff),
									flagOutliers=integer(ly), as.integer(B.init),
									as.double(tol.init), as.integer(trans),
									logLike=as.double(0),
									as.integer(control$B.lambda), as.integer(control$B.brent),
									as.double(control$tol.brent), as.double(control$xLow),
									as.double(control$xUp),
									as.double(t(Mu0)),
									as.double(kappa0),
									as.double(nu0),
									as.double(Lambda0),
									as.double(Omega0),
									as.integer(.model),
									PACKAGE="flowClust"))
				}
				if (class(obj)!="try-error")
				{
					maxLabel[[j]] <- label
					maxLogLike[j] <- obj$logLike
				}
			}
			ind <- order(maxLogLike, decreasing=T, na.last=NA)
		}
	}
	else
	{
		if(usePrior=="yes")
		{
			M<-0
                        #FIXME priors with transformation is currently borked
			#Assuming Mu0 is on transformed scale
			if(lambda!=1){
				ytmp<-box(y,lambda);
			}else{
				ytmp<-y;
			}
      # We use the prior densities to initialize the cluster labels to the
      # mixture component (cluster) that maximizes the prior probability.
      prob <- sapply(seq_len(K), function(k) {
        with(prior, w0[k] * dmvt(x = ytmp, mu = Mu0[k, ], sigma = Lambda0[k,,],
                                 nu = nu0[k], lambda = lambda)$value)
      })
      label <- apply(prob, 1, which.max)
		}
	}


# long EMs
	for (M in ind)
	{
		if (nu != Inf)
		{
			obj <- try(.C("flowClust", as.double(t(y)), as.integer(ly),
							as.integer(py), as.integer(K[i]),
							w=rep(0,K[i]), mu=rep(0,K[i]*py),
							precision=initprec,
							lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
							nu=as.double(rep(nu,K[i])),
							z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
							as.integer(if (M==0) label else maxLabel[[M]]), uncertainty=double(ly),
							as.double(rep(u.cutoff,K[i])), as.double(z.cutoff),
							flagOutliers=integer(ly), as.integer(B),
							as.double(tol), as.integer(trans),
							as.integer(nu.est), logLike=as.double(0),
							as.integer(control$B.lambda), as.integer(control$B.brent),
							as.double(control$tol.brent), as.double(control$xLow),
							as.double(control$xUp), as.double(control$nuLow),
							as.double(control$nuUp),
							mu0=as.double(t(Mu0)),
							as.double(kappa0),
							nu0=as.double(nu0),
							lambda0=as.double(Lambda0),
							omega0=as.double(Omega0),
							w0=as.double(w0),
							as.integer(.model),
							oorder=as.integer(oorder),
							PACKAGE="flowClust"))
							if (class(obj)=="try-error"){
								message("flowClust failed")
							}

		}
		else
		{
			obj <- try(.C("flowClustGaussian", as.double(t(y)), as.integer(ly),
							as.integer(py), as.integer(K[i]),
							w=rep(0,K[i]), mu=rep(0,K[i]*py),
							precision=initprec,
							lambda=as.double(rep(lambda, length.out=(if (trans>1) K[i] else 1))),
							z=rep(0,ly*K[i]), u=rep(0,ly*K[i]),
							as.integer(if (M==0) label else maxLabel[[M]]), uncertainty=double(ly),
							as.double(q.cutoff), as.double(z.cutoff),
							flagOutliers=integer(ly), as.integer(B),
							as.double(tol), as.integer(trans),
							logLike=as.double(0),
							as.integer(control$B.lambda), as.integer(control$B.brent),
							as.double(control$tol.brent), as.double(control$xLow),
							as.double(control$xUp),
							as.double(t(Mu0)),
							as.double(kappa0),
							as.double(nu0),
							as.double(Lambda0),
							as.double(Omega0),
							as.integer(.model),
							PACKAGE="flowClust"))
			if (class(obj)!="try-error"){

				obj$nu <- Inf
			}
		}
		if (class(obj)!="try-error")
			break
	}
	if (class(obj)=="try-error")
		stop(geterrmessage())

	if(usePrior=="vague"){
		trans<-1;
	}
# output obj$precision to sigma
	sigma <- array(0, c(K[i], py, py))
	precision <- matrix(obj$precision, K[i], py * py, byrow=TRUE)
	for (k in 1:K[i])
		sigma[k,,] <- matrix(precision[k,], py, py, byrow = TRUE)

# output BIC & ICL
	BIC <- 2*obj$logLike - log(ly) * (K[i]*(py+1)*py/2 + K[i]*py + K[i]-1 + (if (trans>1) K[i] else trans) + (if (nu.est>1) K[i] else abs(nu.est)))
	z <- matrix(obj$z, ly, K[i], byrow = TRUE)
	ICL <- BIC + 2 * sum(z*log(z), na.rm = TRUE)
  if(is.null(include)){
    # output z, u, label, uncertainty, flagOutliers
    z <- u <- matrix(NA, ly, K[i])
    z <- matrix(obj$z, ly, K[i], byrow=TRUE)
    u <- matrix(obj$u, ly, K[i], byrow=TRUE)
    #cat(M);
    tempLabel <- if (M==0) label else maxLabel[[M]]
    label <- uncertainty <- flagOutliers <- rep(NA, ly)
    label <- tempLabel
    uncertainty <- obj$uncertainty
    flagOutliers <- as.logical(obj$flagOutliers)
  }else{
    # output z, u, label, uncertainty, flagOutliers
    z <- u <- matrix(NA, length(include), K[i])
    z[include,] <- matrix(obj$z, ly, K[i], byrow=TRUE)
    u[include,] <- matrix(obj$u, ly, K[i], byrow=TRUE)
    #cat(M);
    tempLabel <- if (M==0) label else maxLabel[[M]]
    label <- uncertainty <- flagOutliers <- rep(NA, length(include))
    label[include] <- tempLabel
    uncertainty[include] <- obj$uncertainty
    flagOutliers[include] <- as.logical(obj$flagOutliers)
    
  }

# output reordered prior
	prior$Mu0<-matrix({if(all(!is.null(obj$mu0))){obj$mu0}else{NA}},K[i],py,byrow=TRUE);
	prior$Lambda0<-aperm(array({if(all(!is.null(obj$lambda0))){obj$lambda0}else{NA}},c(py,py,K[i])),c(3,1:2))
	prior$Omega0<-aperm(array({if(all(!is.null(obj$omega0))){obj$omega0}else{NA}},c(py,py,K[i])),c(3,1:2))
	prior$nu0<-{if(all(!is.null(obj$nu0))){obj$nu0}else{NA}}
	prior$w0<-{if(all(!is.null(obj$w0))){obj$nu0}else{NA}}
#omit
#result<- new("flowClust", expName=expName, varNames=varNames, K=K[i],
#		w=obj$w, mu=matrix(obj$mu, K[i], py, byrow=TRUE), sigma=sigma,
#		lambda=(if (trans>0) obj$lambda else numeric(0)), nu=(if (nu.est>1) obj$nu else obj$nu[1]), z=z,
#		u=u, label=label, uncertainty=uncertainty,
#		ruleOutliers=ruleOutliers, flagOutliers=flagOutliers, rm.min=sum(rm.min),
#		rm.max=sum(rm.max), logLike=obj$logLike, BIC=BIC, ICL=ICL);
class(prior)<-"list";
prior$order<-obj$oorder;
	if(trans==1&obj$lambda==1){
		obj$mu<-rbox(obj$mu,obj$lambda)
	}
	#do nothing in particular if trans>1
        #Not sure if the above does the right thing when trans=1, and the returned lambda=1, and usePrior="no"
	result<- new("flowClust", expName=expName, varNames=varNames, K=K[i],
			w=obj$w, mu=matrix(obj$mu, K[i], py, byrow=TRUE), sigma=sigma,
			lambda= obj$lambda, nu=(if (nu.est>1) obj$nu else obj$nu[1]), z=z,
			u=u, label=label, uncertainty=uncertainty,
			ruleOutliers=ruleOutliers, flagOutliers=flagOutliers, rm.min=sum(rm.min),
			rm.max=sum(rm.max), logLike=obj$logLike, BIC=BIC, ICL=ICL,prior=prior);
	# if(!any(is.na(prior))&usePrior=="yes"&ruleOutliers[1]==0){
		# label<-.fcbMap(result,ruleOutliers[2])
		# result@flagOutliers<-label==0
		# label[result@flagOutliers]<-NA;
		# result@label<-label;
	# }
	result
}
