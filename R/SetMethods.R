## show and summary methods
setMethod("show", "flowClust",
          function(object)
      {
          cat("Object of class 'flowClust'","\n")
          cat("This object has the following slots: \n")
          cat("expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label,",
              "uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max,",
              "logLike, BIC, ICL\n")
      })


setMethod("summary", "flowClust",
          function(object)
      {
          cat("** Experiment Information ** \n")
          cat("Experiment name:",object@expName,"\n")
          cat("Variables used:",object@varNames,"\n")	
          cat("** Clustering Summary ** \n")	
          cat("Number of clusters:",object@K,"\n")
          cat("Proportions:",object@w,"\n")	
          cat("** Transformation Parameter ** \n")
          cat("lambda:",object@lambda,"\n")  
          cat("** Information Criteria ** \n")
          cat("Log likelihood:",object@logLike,"\n")
          cat("BIC:",object@BIC,"\n")
          cat("ICL:",object@ICL,"\n")
          cat("** Data Quality ** \n")
          if (!is.na(object@rm.max))
              cat("Number of points filtered from above: ",
                  object@rm.max, " (",
                  round(object@rm.max/nrow(object@z)*100, 2), "%)\n", sep="")
          if (!is.na(object@rm.min))
              cat("Number of points filtered from below: ", object@rm.min,
                  " (", round(object@rm.min/nrow(object@z)*100, 2), "%)\n",
                  sep="")
          ruleOutliers(object)
          cat("Number of outliers: ", sum(object@flagOutliers, na.rm=T), " (",
              round(sum(object@flagOutliers, na.rm=T)/nrow(object@z)*100, 2),
              "%)\n", sep="")
          cat("Uncertainty summary: \n")
          summary(object@uncertainty)
      })


setMethod("show", signature("tmixFilter"),
          function(object)
      {
          cat("A t-mixture filter named '", object@filterId, "'", sep="")
          if (length(parameters(object))==1 && parameters(object)=="")
              cat("\n")
          else
              cat(":", paste(parameters(object), collapse=", "), "\n")
          cat("** Clustering Settings ** \n")	
          cat("Number of clusters:",object@K,"\n")
          cat("Degrees of freedom of t distribution:", object@nu, "\n")
          cat("Maximum number of EM iterations:", object@B, "\n")
          if (object@trans) {
              cat("Transformation selection will be performed")
              if (object@lambda!=1)
                  cat(" (Initial transformation: lambda =", object@lambda,
                      "\n")
              else
                  cat("\n")
          }
          cat("\n")
          cat("Rule of identifying outliers: ")
          if (is.na(object@u.cutoff))
              cat(round(object@level*100, 2), "% quantile", sep="")
          else
              cat("Cutoff at", round(object@u.cutoff, 2))
          if (object@z.cutoff==0)
              cat("\n")
          else
              cat(", probability of assignment <", round(object@z.cutoff, 2),
                  "\n")
          })


setMethod("show", signature("tmixFilterResult"),
          function(object)
      {
          cat("Object of class 'tmixFilterResult'\n")
          cat("This object stores the result of applying a t-mixture filter",
              "named '", object@filterId, "' on ", object@frameId,
              ", which has an experiment name '", object@expName, "'\n", sep="")
              cat("The following parameters have been used in filtering:\n")
              cat("\t", paste(object@varNames, collapse=", "), "\n", sep="")
              cat("\n")
              cat("The 'tmixFilterResult' class extends the 'flowClust'",
                  "class, which contains the following slots:\n")
          cat("expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label,",
              "uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max,",
              "logLike, BIC, ICL\n")
          })


setMethod("summary", signature("tmixFilterResult"),
          function(object)
      {
          selectMethod("summary", "flowClust")(object)
      })






## %in% methods
## This one is quite weird if you think of the semantics it implies...
setMethod("%in%", signature("ANY", "flowClust"),
          function(x, table)
      {
          !table@flagOutliers & !is.na(table@flagOutliers)
      })


setMethod("%in%", signature("flowFrame", "tmixFilterResult"),
          function(x, table)
      {
          selectMethod("%in%", c("ANY", "flowClust"))(x, table)
      })


setMethod("%in%", signature("flowFrame", "tmixFilter"),
          function(x, table)       
      {
          vn <- if (table@parameters[1]=="") NULL else table@parameters
          min <- if (is.na(table@min[1])) NULL else table@min
          max <- if (is.na(table@max[1])) NULL else table@max
          uc <- if (is.na(table@u.cutoff)) NULL else table@u.cutoff
          flowClust(x, expName=table@expName, varNames=vn,
                    K=table@K, B=table@B, tol=table@tol,
                    nu=table@nu, lambda=table@lambda,
                    trans=table@trans, min.count=table@min.count,
                    max.count=table@max.count, min=min, max=max,
                    level=table@level, u.cutoff=uc,
                    z.cutoff=table@z.cutoff,
                    randomStart=table@randomStart, seed=table@seed)
      })


setMethod("%in%", signature("ANY", "tmixFilterResult"),
          function(x, table)
      {
          callNextMethod()
      })






## Object coercion
setAs("flowClust","logical", function(from) 1 %in% from )


setAs("tmixFilterResult","logical", function(from) 1 %in% from)


setAs("flowClust","filterResult", function(from)
      new("tmixFilterResult", from))


setAs("flowClust","tmixFilterResult", function(from)
      new("tmixFilterResult", from, subSet=factor(Map(from))))



## Subsetting and splitting
setMethod("[",
          signature=signature("flowFrame","flowClust"),
          definition=function(x,i,j,...,drop=FALSE)
      {
          i <- as(i, "tmixFilterResult")
          callGeneric(x,i,j,...,drop=drop)
      })


setMethod("[",
          signature=signature("flowFrame","tmixFilterResult"),
          definition=function(x, i, j, ..., drop=FALSE)
      {
          selectMethod("[",
                       c("flowFrame", "filterResult"))(x, i, j, ..., drop=drop)
      })


setMethod("Subset", signature("flowFrame","flowClust"),
          function(x,subset,...)
      {
          subset <- as(subset, "tmixFilterResult")
          callGeneric()
      })

setMethod("Subset",
          signature("flowFrame", "tmixFilterResult"),
          function(x,subset,...)
      {
          selectMethod("Subset",
                       c("flowFrame", "filterResult"))(x,subset,...)
      })


setMethod("Subset", signature("data.frame", "flowClust"),
          function(x, subset, ..., select)
      {
          object <- 
              if (missing(select))
                  as.data.frame(x[x %in% subset, , drop=FALSE])
              else
                  as.data.frame(x[x %in% subset, select, drop=FALSE])
          return(object)
      })


setMethod("Subset", signature("matrix","flowClust"),
          function(x, subset, ..., select)
      {
          x <- as.data.frame(x)
          callGeneric()
      })


setMethod("Subset", signature("vector","flowClust"),
          function(x,subset,...)
      {
          object <- x[x %in% subset]       
          return(object)
      })



.spltVsPop <- function(pop, splt, f){
    DEPR_MESSAGE <- paste("The 'split' argument is deprecated.\nPlease use",
                          "'population' instead.")
    if(is.null(pop)){
        if(!is.null(splt)){
            pop <- splt
            warning(msg=DEPR_MESSAGE, call.=FALSE)
        }else{
            pop <- as.list(1:f@K)
        }
    }else if(!is.null(splt)){
        warning("Both arguments 'population' and 'split' are specified.",
                "\n'split' is deprecated and will be ignored.", call.=FALSE)
    }
    return(pop)
}


setMethod("split",
          signature(x="data.frame", f="flowClust", drop="ANY"),
          function(x, f, drop=FALSE, population=NULL, split=NULL,
                   rm.outliers=TRUE)
      {
          population <- .spltVsPop(population, split, f)
          object <- vector("list", length(population))
          for (i in 1:length(population)) 
              object[[i]] <- x[is.element(Map(f, rm.outliers),
                                          population[[i]]),, drop=FALSE]
          names(object) <- names(population)
          return(object)
      })  


setMethod("split",
          signature(x="matrix", f="flowClust", drop="ANY"),
          function(x, f, drop=FALSE, population=NULL, split=NULL,
                   rm.outliers=TRUE)
      {
          split(as.data.frame(x), f, drop=drop, population=population,
                split=split, rm.outliers=rm.outliers)
      })


setMethod("split",
          signature(x="vector", f="flowClust", drop="ANY"),
          function(x, f, drop=FALSE, population=NULL, split=NULL,
                   rm.outliers=TRUE)
      {
          population <- .spltVsPop(population, split, f)
          object <- vector("list", length(population))
          for (i in 1:length(population))
              object[[i]] <- x[is.element(Map(f, rm.outliers), population[[i]])]
          names(object) <- names(population)
          return(object)
      })


setMethod("split",
          signature(x="flowFrame", f="flowClust", drop="ANY"),
          function(x, f, drop=FALSE, population=NULL, split=NULL,
                   rm.outliers=TRUE, ...)
      {
          f <- as(f, "tmixFilterResult")
          selectMethod("split", c("flowFrame", "tmixFilterResult"))(x,f,drop,population,split,rm.outliers, ...)
      })


setMethod("split",
          signature(x="flowFrame", f="tmixFilterResult", drop="ANY"),
          function(x, f, drop=FALSE, population=NULL, split=NULL,
                   rm.outliers=TRUE, ...)
      {
          population <- sapply(.spltVsPop(population, split, f), as.character)
          if(!is.list(population))
              population <- as.list(population)
          subSet <- factor(Map(f, rm.outliers))
          x <- x[!is.na(subSet),]
          subSet <- subSet[!is.na(subSet)]
          f@subSet <- subSet
          selectMethod("split",
                     c("flowFrame", "multipleFilterResult"))(x, f, drop, population=population, ...)
      })



## Misc methods
setGeneric("ruleOutliers", function(object) standardGeneric("ruleOutliers"))


setMethod("ruleOutliers",
          signature("flowClust"),
          function(object)
      {
          cat("Rule of identifying outliers: ")
          if (object@ruleOutliers[1]==0)
              cat(round(object@ruleOutliers[2]*100, 2), "% quantile", sep="")
          else
              cat("Cutoff at", round(object@ruleOutliers[2], 2))
          if (object@ruleOutliers[3]==0)
              cat("\n")
          else
              cat(",\n                              probability of ",
                  "assignment <", round(object@ruleOutliers[3], 2), "\n")
      })


setGeneric("ruleOutliers<-",
           function(object,value) standardGeneric("ruleOutliers<-"))


setReplaceMethod("ruleOutliers",
                 signature("flowClust","list"),
                 function(object,value=list(level=NULL, u.cutoff=NULL,
                                 z.cutoff=NULL))
             {
                 if (is.null(value$u.cutoff) && !is.null(value$level))
                     object@ruleOutliers[1:2] <- c(0, value$level)
                 else if (!is.null(value$u.cutoff))
                     object@ruleOutliers[1:2] <- c(1, value$u.cutoff)
                 if (!is.null(value$z.cutoff))
                     object@ruleOutliers[3] <- value$z.cutoff
                 if (object@nu!=Inf) {
                     if (object@ruleOutliers[1]==0){     # 0 means quantile
                         py <- ncol(object@mu)
                         cc <- py * qf(object@ruleOutliers[2], py, object@nu)
                         u.cutoff <- (object@nu+py) / (object@nu+cc)
                         result <- (importance(object, assign=T) < u.cutoff)
                     }else{
                         result <- (importance(object, assign=T) <
                                    object@ruleOutliers[2])
                     }
                 } else {
                     py <- ncol(object@mu)
                     q.cutoff <- qchisq(object@ruleOutliers[2], py)
                     result <- ( importance(object, assign=T) > q.cutoff)
                 }
                 if (object@ruleOutliers[3]>0)
                     result <- result | (posterior(object, assign=T) <
                                         object@ruleOutliers[3])
                 object@flagOutliers <- result
                 ruleOutliers(object)
                 object
             })


setGeneric("Map")


setMethod("Map",
          signature(f="flowClust"),
          function(f, rm.outliers=TRUE, ...)
      {
          result <- map(f@z)
          if (rm.outliers) result[which(f@flagOutliers)] <- NA
          result
      })






## helper functions

criterion <- function(object, type="BIC")
{
    if (type=="BIC") object@BIC  else if (type=="ICL") object@ICL  else if (type=="logLike") object@logLike
}


posterior <- function(object, assign=FALSE)
{
    if (!assign) object@z  else {
        result <- rep(NA, nrow(object@z))
        result[!is.na(object@flagOutliers)] <-
            t(object@z)[t(unmap(map(object@z))==T)]
        result
    }
}


importance <- function(object, assign=FALSE)
{
    if (!assign) object@u  else {
        result <- rep(NA, nrow(object@u))
        result[!is.na(object@flagOutliers)] <-
            t(object@u)[t(unmap(map(object@z))==T)]
        result
    }
}


uncertainty <- function(object)  object@uncertainty


getEstimates <- function(object, data)
{
    if (object@lambda==1)
        list(proportions=object@w, locations=object@mu,
             dispersion=object@sigma)
    else{
        if (missing(data))
            list(proportions=object@w, locations=rbox(object@mu,
                                       object@lambda))
        else
        {
            if(is(data,"flowFrame"))
            {
                y<-as.matrix(exprs(data)[,object@varNames])
            }
            else if(is(data,"matrix"))
            {
                if(object@varNames=="Not Available")
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
            if (object@nu!=Inf){
                obj <- .C("getEstimates", as.double(t(y)), as.integer(ly),
                          as.integer(py), as.integer(K), mu=rep(0.0,K*py),
                          precision=rep(0.0,K*py*py), as.double(object@nu),
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
            list(proportions=object@w, locations=rbox(object@mu,
                                       object@lambda),
                 locationsC=matrix(obj$mu,K,py,byrow=TRUE), dispersion=sigma)
        }
    }
}




setMethod("plot",
          signature("flowFrame", "tmixFilterResult"),
          function(x, y, z=NULL, ...)
      {
          if(is.null(z))
              z <- y@varNames
          if(length(z)>2) {
              warning("Only plotting the first two parameters.")
              z <- z[1:2]
          }
          if (length(z)==2)
              selectMethod("plot",
                           signature("flowClust", "missing"))(x=y, data=x, subset=z, ...)
          else
              selectMethod("hist",
                           signature("flowClust"))(x=y, data=x, subset=z, ...)
    }
)
