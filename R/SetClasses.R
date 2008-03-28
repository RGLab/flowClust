setClass("flowClust",
         representation(expName="character",varNames="character",K="numeric",mu="matrix", sigma="array",w="vector",z="matrix",u="matrix",label="vector",uncertainty="vector",ruleOutliers="vector",flagOutliers="vector",rm.max="numeric",rm.min="numeric",lambda="numeric",nu="numeric",logLike="numeric",BIC="numeric",ICL="numeric"),
         prototype(expName=character(0),varNames=character(0),K=numeric(0),mu=matrix(numeric(0), nrow=0, ncol=0),sigma=array(numeric(0), c(0,0,0)),w=rep(numeric(0),0),
		 z=matrix(numeric(0), nrow=0, ncol=0),u=matrix(numeric(0), nrow=0, ncol=0),
		label=matrix(numeric(0), nrow=0, ncol=0),
		uncertainty=rep(numeric(0), 0), ruleOutliers=rep(numeric(0),0),
		flagOutliers=rep(logical(0), 0), rm.max=numeric(0), rm.min=numeric(0),
		lambda=double(1),nu=numeric(0),logLike=numeric(0),BIC=numeric(0),ICL=numeric(0)))
		
setClass("flowDens", representation(dx="matrix", dy="matrix", value="matrix"), prototype(dx=matrix(numeric(0), nrow=0, ncol=1), dy=matrix(numeric(0), nrow=0, ncol=1), value=matrix(numeric(0), nrow=0, ncol=0)))


setClass("tmixFilter",
    representation(expName="character", K="numeric", B="numeric", tol="numeric", nu="numeric", lambda="numeric", trans="logical", min.count="numeric", max.count="numeric", min="vector", max="vector", level="numeric", u.cutoff="numeric", z.cutoff="numeric", randomStart="logical"),
    prototype(expName="Flow Experiment", K=numeric(0), B=500, tol=1e-5, nu=4, lambda=1, trans=TRUE, min.count=10, max.count=10, min=NA, max=NA, level=0.9, u.cutoff=NA_real_, z.cutoff=0, randomStart=FALSE),
    contains="parameterFilter")

tmixFilter <- function(filterId="tmixFilter", parameters="", ...)
{
    if (!"K" %in% names(list(...))) stop("The number of clusters, K, must be provided.")
    new("tmixFilter", filterId=filterId, parameters=parameters, ...)
}

setClass("tmixFilterResult", contains=c("flowClust", "filterResult"))

setMethod("filter", signature(x="flowFrame", filter="tmixFilter"),
    function(x, filter) {
        result <- flowClust(x, expName=filter@expName, varNames=(if (filter@parameters[1]=="") NULL else filter@parameters), K=filter@K, B=filter@B, tol=filter@tol, nu=filter@nu, lambda=filter@lambda, trans=filter@trans, min.count=filter@min.count, max.count=filter@max.count, min=(if (is.na(filter@min[1])) NULL else filter@min), max=(if (is.na(filter@max[1])) NULL else filter@max), level=filter@level, u.cutoff=(if (is.na(filter@u.cutoff)) NULL else filter@u.cutoff), z.cutoff=filter@z.cutoff, randomStart=filter@randomStart)
        result <- new("tmixFilterResult", result, filterId=identifier(filter), frameId=identifier(x))
        filterDetails(result, identifier(filter)) <- filter
        result
    }
)


setMethod("plot", signature("flowFrame", "tmixFilterResult"), function(x, y, z=NULL, ...) {
    if(is.null(z))	z = y@varNames
    if(length(z)>2) {
        warning("Only plotting the first two parameters.")
        z = z[1:2]
    }
    if (length(z)==2) selectMethod("plot", signature("flowClust", "missing"))(x=y, data=x, subset=z, ...)  else selectMethod("hist", signature("flowClust"))(x=y, data=x, subset=z, ...)
    }
)
