setClass("flowClust",
         representation(expName="character", varNames="character",
                        K="numeric", w="vector", mu="matrix",
                        sigma="array", lambda="numeric", nu="numeric",
                        z="matrix", u="matrix", label="vector",
                        uncertainty="vector", ruleOutliers="vector",
                        flagOutliers="vector", rm.min="numeric",
                        rm.max="numeric", logLike="numeric", BIC="numeric",
                        ICL="numeric"),
         prototype(expName=character(0), varNames=character(0), K=numeric(0),
                   w=rep(numeric(0),0), mu=matrix(numeric(0), nrow=0, ncol=0),
                   sigma=array(numeric(0), c(0,0,0)), lambda=numeric(0),
                   nu=numeric(0), z=matrix(numeric(0), nrow=0, ncol=0),
                   u=matrix(numeric(0), nrow=0, ncol=0),
                   label=rep(numeric(0),0), uncertainty=rep(numeric(0),0),
                   ruleOutliers=rep(numeric(0),0),
                   flagOutliers=rep(logical(0), 0), rm.min=numeric(0),
                   rm.max=numeric(0), logLike=numeric(0), BIC=numeric(0),
                   ICL=numeric(0)))


setClass("flowClustList",
         representation("list", index="numeric", criterion="character"),
         prototype(vector("list",0), index=numeric(0), criterion="BIC"),
         validity=function(object)
                  {
                      if (!all(sapply(object@.Data, is, class2="flowClust")))
                          return("Not a list of flowClust results!")
                      else return(TRUE)
                  })


setClass("flowDens",
         representation(dx="matrix", dy="matrix", value="matrix"),
         prototype(dx=matrix(numeric(0), nrow=0, ncol=1),
                   dy=matrix(numeric(0), nrow=0, ncol=1),
                   value=matrix(numeric(0), nrow=0, ncol=0)))


setClass("tmixFilter",
         representation(expName="character", K="numeric", B="numeric",
                        tol="numeric", nu="numeric", lambda="numeric",
                        nu.est="numeric", trans="numeric", min.count="numeric",
                        max.count="numeric", min="vector", max="vector",
                        level="numeric", u.cutoff="numeric", z.cutoff="numeric",
                        randomStart="numeric", B.init="numeric", tol.init="numeric",
                        seed="numeric", criterion="character", control="list"),
         prototype(expName="Flow Experiment", K=numeric(0), B=500, tol=1e-5,
                   nu=4, lambda=1, nu.est=0, trans=1, min.count=10, max.count=10,
                   min=NA, max=NA, level=0.9, u.cutoff=NA_real_, z.cutoff=0,
                   randomStart=10, B.init=500, tol.init=1e-2, seed=1, criterion="BIC", 
                   control=vector("list",0)),
         contains="parameterFilter")


tmixFilter <- function(filterId="tmixFilter", parameters="", ...)
{
    if (!"K" %in% names(list(...)))
        stop("The number of clusters, K, must be provided.")
    new("tmixFilter", filterId=filterId, parameters=as.list(parameters), ...)
}


setClass("tmixFilterResult", contains=c("flowClust", "multipleFilterResult"))


setClass("tmixFilterResultList",
         contains=c("flowClustList", "multipleFilterResult"))
