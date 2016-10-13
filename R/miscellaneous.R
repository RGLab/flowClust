#' Various Functions for Retrieving Information from Clustering Results
#' 
#' Various functions are available to retrieve the information criteria
#' (\code{criterion}), the posterior probabilities of clustering memberships
#' \eqn{z} (\code{posterior}), the \dQuote{weights} \eqn{u}
#' (\code{importance}), the uncertainty (\code{uncertainty}), and the estimates
#' of the cluster proportions, means and variances (\code{getEstimates})
#' resulted from the clustering (filtering) operation.
#' 
#' These functions are written to retrieve various slots contained in the
#' object returned from the clustering operation.  \code{criterion} is to
#' retrieve \code{object@BIC}, \code{object@ICL} or \code{object@logLike}.  It
#' replacement method modifies \code{object@index} and \code{object@criterion}
#' to select the best model according to the desired criterion.
#' \code{posterior} and \code{importance} provide a means to conveniently
#' retrieve information stored in \code{object@z} and \code{object@u}
#' respectively.  \code{uncertainty} is to retrieve \code{object@uncertainty}.
#' \code{getEstimates} is to retrieve information stored in \code{object@mu}
#' (transformed back to the original scale) and \code{object@w}; when the data
#' object is provided, an approximate variance estimate (on the original scale,
#' obtained by performing one M-step of the EM algorithm without taking the
#' Box-Cox transformation) will also be computed.
#' 
#' @aliases criterion criterion,flowClust-method criterion,flowClustList-method
#' criterion<- criterion<-,flowClustList,character-method posterior importance
#' uncertainty getEstimates
#' @param object Object returned from \code{\link{flowClust}} or
#' \code{\link[=tmixFilter]{filter}}.  For the replacement method of
#' \code{criterion}, the object must be of class \code{\link{flowClustList}} or
#' \code{\link{tmixFilterResultList}}.
#' @param \dots Further arguments. Currently this is \code{type}, a character
#' string.  May take \code{"BIC"}, \code{"ICL"} or \code{"logLike"}, to specify
#' the criterion desired.
#' @param type,value A character string stating the criterion used to choose the
#' best model.  May take either \code{"BIC"} or \code{"ICL"}.
#' @param max whether \code{criterion} should return the max value
#' @param show.K whether \code{criterion} should return K
#' @param assign A logical value.  If \code{TRUE}, only the quantity (\code{z}
#' for \code{posterior} or \code{u} for \code{importance}) associated with the
#' cluster to which an observation is assigned will be returned.  Default is
#' \code{FALSE}, meaning that the quantities associated with all the clusters
#' will be returned.
#' @param data A numeric vector, matrix, data frame of observations, or object
#' of class \code{flowFrame}; an optional argument.  This is the object on
#' which \code{\link{flowClust}} or \code{\link[=tmixFilter]{filter}} was
#' performed.
#' @return Denote by \eqn{K} the number of clusters, \eqn{N} the number of
#' observations, and \eqn{P} the number of variables.  For \code{posterior} and
#' \code{importance}, a matrix of size \eqn{N \times K}{N x K} is returned if
#' \code{assign=FALSE} (default).  Otherwise, a vector of size \eqn{N} is
#' outputted.  \code{uncertainty} always outputs a vector of size \eqn{N}.
#' \code{getEstimates} returns a list with named elements, \code{proportions},
#' \code{locations} and, if the data object is provided, \code{dispersion}.
#' \code{proportions} is a vector of size \eqn{P} and contains the estimates of
#' the \eqn{K} cluster proportions.  \code{locations} is a matrix of size
#' \eqn{K \times P}{K x P} and contains the estimates of the \eqn{K} mean
#' vectors transformed back to the original scale (i.e., \code{rbox(object@mu,
#' object@lambda)}).  \code{dispersion} is an array of dimensions \eqn{K \times
#' P \times P}{K x P x P}, containing the approximate estimates of the \eqn{K}
#' covariance matrices on the original scale.
#' @note When \code{object@nu=Inf}, the Mahalanobis distances instead of the
#' \dQuote{weights} are stored in \code{object@u}.  Hence, \code{importance}
#' will retrieve information corresponding to the Mahalanobis distances. % If
#' the \code{assign} argument is set to \code{TRUE}, only the quantities
#' corresponding to assigned observations will be returned.  Quantities
#' corresponding to unassigned observations (outliers and filtered
#' observations) will be reported as \code{NA}.  Hence, A change in the rule to
#' call outliers will incur a change in the number of \code{NA} values
#' returned.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}},
#' \code{\link{Map}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords cluster
#' @rdname miscellaneous
#' @export 
setGeneric("criterion", function(object, ...) standardGeneric("criterion"))
#' @rdname miscellaneous
setMethod("criterion", signature(object="flowClust"),
    function(object, type="BIC")
    {
      if (type=="BIC") object@BIC  else if (type=="ICL") object@ICL  else if (type=="logLike") object@logLike
    })
#' @rdname miscellaneous
setMethod("criterion", signature(object="flowClustList"),
    function(object, type="BIC", max=FALSE, show.K=FALSE)
    {
      values <- rep(0, length(object))
      if (show.K) K <- rep(0, length(object))
      for (i in 1:length(object))
      {
        if (type=="BIC") values[i] <- object[[i]]@BIC else
        if (type=="ICL") values[i] <- object[[i]]@ICL else
        if (type=="logLike") values[i] <- object[[i]]@logLike
        if (show.K) K[i] <- object[[i]]@K     
      }
      if (max)
      {
        if (show.K) K <- K[which.max(values)] else values <- max(values)
      }
      if (show.K) K else values
    })

#' @rdname miscellaneous
#' @export 
setGeneric("criterion<-", function(object,value) standardGeneric("criterion<-"))
#' @rdname miscellaneous
setReplaceMethod("criterion", signature("flowClustList","character"),
    function(object,value)
    {
      values <- criterion(object, value)
      object@index <- which.max(values)
      object@criterion <- value
      object
    })


#' @rdname miscellaneous
#' @export 
posterior <- function(object, assign=FALSE)
{
  if (is(object,"flowClustList")) object <- object[[object@index]]
  if (!assign) object@z  else {
    result <- rep(NA, nrow(object@z))
    result[!is.na(object@flagOutliers)] <-
        apply(object@z, 1, max)[!is.na(object@flagOutliers)]
#            t(object@z)[t(unmap(map(object@z))==T)]
    result
  }
}

#' @rdname miscellaneous
#' @export 
importance <- function(object, assign=FALSE)
{
  if (is(object,"flowClustList")) object <- object[[object@index]]
  if (!assign) object@u  else {
    result <- rep(NA, nrow(object@u))
    result[!is.na(object@flagOutliers)] <-
        object@u[cbind(1:nrow(object@u), max.col(object@z, "first"))][!is.na(object@flagOutliers)]
#            t(object@u)[t(unmap(map(object@z))==T)]
    result
  }
}

#' @rdname miscellaneous
#' @export 
uncertainty <- function(object)
{
  if (is(object,"flowClustList")) object <- object[[object@index]]
  object@uncertainty
}

#' @rdname miscellaneous
#' @export 
getEstimates <- function(object, data)
{
  if (is(object,"flowClustList")) object <- object[[object@index]]
  if (length(object@lambda)==0)
    list(proportions=object@w, locations=object@mu,
        dispersion=object@sigma)
  else{
    if (missing(data))
      list(proportions=object@w, locations={if((object@lambda!=1))
            {
              rbox(object@mu,object@lambda)
            }else{
              object@mu
            }
          })
    
    else
    {
      if(is(data,"flowFrame"))
      {
        y<-as.matrix(exprs(data)[,object@varNames])
      }
      else if(is(data,"matrix"))
      {
        if(object@varNames[1]=="Not Available")
          y<-data
        else
          y<-as.matrix(data[,object@varNames])
      }
      else if(is(data,"data.frame"))
      {
        y<-as.matrix(data[,object@varNames])
      }
      else if(is(data,"vector"))
      {
        y<-matrix(data)
      }    
      include <- !is.na(object@uncertainty)
      y <- as.matrix(y[include,])
      z <- as.matrix(object@z[include,])
      u <- as.matrix(object@u[include,])
      ly<-nrow(y)
      py<-ncol(y)
      K<-object@K
      if (all(object@nu!=Inf)){
        obj <- .C("getEstimates", as.double(t(y)), as.integer(ly),
            as.integer(py), as.integer(K), mu=rep(0.0,K*py),
            precision=rep(0.0,K*py*py), as.double(rep(object@nu,length.out=K)),
            as.double(t(z)), as.double(t(u)))
      }else{
        obj <- .C("getEstimatesGaussian", as.double(t(y)),
            as.integer(ly), as.integer(py), as.integer(K),
            mu=rep(0.0,K*py), precision=rep(0.0,K*py*py),
            as.double(t(z)))
      }
      sigma <- array(0,c(K,py,py))
      precision <- matrix(obj$precision, K, py*py, byrow=TRUE)
      for(k in 1:K) sigma[k,,] <- matrix(precision[k,], py, py,
            byrow=TRUE)
      list(proportions=object@w, {if((object@lambda!=1)){
              locations=rbox(object@mu,
                  object@lambda)
            }else{
              object@mu
            }},
          locationsC=matrix(obj$mu,K,py,byrow=TRUE), dispersion=sigma)
    }
  }
}

