setMethod("show", "flowClust",
function(object)
{
    cat("Object of class 'flowClust'","\n")
    cat("This object has the following slots: \n")
    cat("expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label, uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max, logLike, BIC, ICL\n")
}
)

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
    if (!is.na(object@rm.max)) cat("Number of points filtered from above: ", object@rm.max, " (", round(object@rm.max/nrow(object@z)*100, 2), "%)\n", sep="")
    if (!is.na(object@rm.min)) cat("Number of points filtered from below: ", object@rm.min, " (", round(object@rm.min/nrow(object@z)*100, 2), "%)\n", sep="")
    ruleOutliers(object)
    cat("Number of outliers: ", sum(object@flagOutliers, na.rm=T)," (",  round(sum(object@flagOutliers, na.rm=T)/nrow(object@z)*100, 2), "%)\n", sep="")
    cat("Uncertainty summary: \n")
    summary(object@uncertainty)
}
)


setMethod("show", signature("tmixFilter"),
    function(object) {
        cat("A t-mixture filter named '", object@filterId, "'", sep="")
        if (length(parameters(object))==1 && parameters(object)=="") cat("\n")
        else cat(":", paste(parameters(object), collapse=", "), "\n")

        cat("** Clustering Settings ** \n")	
        cat("Number of clusters:",object@K,"\n")
        cat("Degrees of freedom of t distribution:", object@nu, "\n")
        cat("Maximum number of EM iterations:", object@B, "\n")
        if (object@trans) {
            cat("Transformation selection will be performed")
            if (object@lambda!=1) cat(" (Initial transformation: lambda =", object@lambda, "\n")  else cat("\n")
        }
        cat("\n")
        cat("Rule of identifying outliers: ")
        if (is.na(object@u.cutoff)) cat(round(object@level*100, 2), "% quantile", sep="")  else cat("Cutoff at", round(object@u.cutoff, 2))
        if (object@z.cutoff==0) cat("\n")  else cat(", probability of assignment <", round(object@z.cutoff, 2), "\n")
    }
)

setMethod("show", signature("tmixFilterResult"),
    function(object) {
        cat("Object of class 'tmixFilterResult'\n")
        cat("This object stores the result of applying a t-mixture filter named '", object@filterId, "' on ", object@frameId, ", which has an experiment name '", object@expName, "'\n", sep="")

        cat("The following parameters have been used in filtering:\n")
        cat("\t", paste(object@varNames, collapse=", "), "\n", sep="")
        cat("\n")
        cat("The 'tmixFilterResult' class extends the 'flowClust' class, which contains the following slots:\n")
        cat("expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label, uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max, logLike, BIC, ICL\n")
    }
)

setMethod("summary", signature("tmixFilterResult"), function(object) selectMethod("summary", "flowClust")(object)
)


setMethod("%in%", signature("ANY", "flowClust"),
    function(x, table){
        !table@flagOutliers & !is.na(table@flagOutliers)
    }
)

setMethod("%in%", signature("flowFrame", "tmixFilterResult"),
    function(x, table){
        selectMethod("%in%", signature("ANY", "flowClust"))(x, table)
    }
)

setMethod("%in%", signature("ANY", "tmixFilterResult"),
    function(x, table){
        selectMethod("%in%", signature("ANY", "flowClust"))(x, table)
    }
)

setAs("flowClust","logical", function(from) 1 %in% from )

setAs("tmixFilterResult","logical", function(from) 1 %in% from)

setMethod("[", signature=signature("flowFrame","flowClust"),
          definition=function(x,i,j,...,drop=FALSE) {
            if(missing(j))
              x[x %in% i,,...,drop=drop]
            else
              x[x %in% i,j,...,drop=drop]
          })

setMethod("[", signature=signature("flowFrame","tmixFilterResult"),
          definition=function(x,i,j,...,drop=FALSE) {
            if(missing(j))
              x[x %in% i,,...,drop=drop]
            else
              x[x %in% i,j,...,drop=drop]
          })


setMethod("Subset", signature("ANY","flowClust"),
    function(x,subset,...){
        if (is(x, "flowFrame")) object <- selectMethod("Subset", signature("flowFrame", "filter"))(x,subset,...)  else
        if (is(x, "matrix")) {
            if (missing(select)) object <- as.matrix(x[x %in% subset,]) else object <- as.matrix(x[x %in% subset, select])
            if (sum(x %in% subset)==1) object <- t(object)
            if (ncol(object)==1) colnames(object) <- if (missing(select)) colnames(x) else select
        }  else
        if (is(x, "data.frame")) {
            if (missing(select)) object <- as.data.frame(x[x %in% subset,]) else object <- as.data.frame(x[x %in% subset, select])
            if (ncol(object)==1) colnames(object) <- if (missing(select)) colnames(x) else select
        }  else
        if (is(x, "vector")) {
            object <- x[x %in% subset]       
        }
        object
    }
)

setMethod("Subset", signature("flowFrame", "tmixFilterResult"),
    function(x,subset,...){
        selectMethod("Subset", signature("ANY", "flowClust"))(x,subset,...)
    }
)


setMethod("split", signature(x="ANY", f="flowClust", drop="missing"),
    function(x, f, split, select, rm.outliers=TRUE) {
        if (missing(split)) split <- as.list(1:f@K)
        object <- vector("list", length(split))
        for (i in 1:length(split)) {
            if (is(x, "flowFrame")) {
                if (missing(select)) object[[i]] <- x[is.element(Map(f, rm.outliers), split[[i]])]  else object[[i]] <- x[is.element(Map(f, rm.outliers), split[[i]]), select]
            }  else
            if (is(x, "matrix")) {
                if (missing(select)) object[[i]] <- as.matrix(x[is.element(Map(f, rm.outliers), split[[i]])])  else object[[i]] <- as.matrix(x[is.element(Map(f, rm.outliers), split[[i]]), select])
                if (sum(is.element(Map(f, rm.outliers), split[[i]]))==1) object[[i]] <- t(object[[i]])
                if (ncol(object[[i]])==1) colnames(object[[i]]) <- if (missing(select)) colnames(x) else select
            }  else
            if (is(x, "data.frame")) {
                if (missing(select)) object[[i]] <- as.data.frame(x[is.element(Map(f, rm.outliers), split[[i]])])  else object[[i]] <- as.data.frame(x[is.element(Map(f, rm.outliers), split[[i]]), select])
                if (ncol(object[[i]])==1) colnames(object[[i]]) <- if (missing(select)) colnames(x) else select
            }
            if (is(x, "vector")) {
                object[[i]] <- x[is.element(Map(f, rm.outliers), split[[i]])]
            }            
        }
        names(object) <- names(split)
        object
    }
)


setMethod("split", signature(x="flowFrame", f="tmixFilterResult", drop="missing"),
    function(x, f, ...) {
        selectMethod("split", signature(x="ANY", f="flowClust", drop="missing"))(x, f, ...)
    }
)


setGeneric("ruleOutliers", function(object) standardGeneric("ruleOutliers"))

setMethod("ruleOutliers", signature("flowClust"),
    function(object) {
        cat("Rule of identifying outliers: ")
        if (object@ruleOutliers[1]==0) cat(round(object@ruleOutliers[2]*100, 2), "% quantile", sep="")  else cat("Cutoff at", round(object@ruleOutliers[2], 2))
        if (object@ruleOutliers[3]==0) cat("\n")  else cat(",\n                              probability of assignment <", round(object@ruleOutliers[3], 2), "\n")
    }
)

setGeneric("ruleOutliers<-", function(object,value) standardGeneric("ruleOutliers<-"))

setReplaceMethod("ruleOutliers", signature("flowClust","list"),
    function(object,value=list(level=NULL, u.cutoff=NULL, z.cutoff=NULL)) {
	    if (is.null(value$u.cutoff) && !is.null(value$level)) object@ruleOutliers[1:2] <- c(0, value$level)  else if (!is.null(value$u.cutoff)) object@ruleOutliers[1:2] <- c(1, value$u.cutoff)
	    if (!is.null(value$z.cutoff)) object@ruleOutliers[3] <- value$z.cutoff

        if (object@nu!=Inf) {
            if (object@ruleOutliers[1]==0) {     # 0 means quantile
                py <- ncol(object@mu)
                cc <- py * qf(object@ruleOutliers[2], py, object@nu)
                u.cutoff <- (object@nu+py) / (object@nu+cc)
                result <- (importance(object, assign=T) < u.cutoff)        
            }  else result <- (importance(object, assign=T) < object@ruleOutliers[2])
        } else {
            py <- ncol(object@mu)
            q.cutoff <- qchisq(object@ruleOutliers[2], py)
            result <- (importance(object, assign=T) > q.cutoff)
        }

        if (object@ruleOutliers[3]>0) result <- result | (posterior(object, assign=T) < object@ruleOutliers[3])
        object@flagOutliers <- result
        ruleOutliers(object)
        object
    }
)


setGeneric("Map")

setMethod("Map", signature(f="flowClust"),
    function(f, rm.outliers=TRUE, ...) {
        result <- map(f@z)
        if (rm.outliers) result[which(f@flagOutliers)] <- NA
        result
    }
)

posterior <- function(object, assign=FALSE) {
    if (!assign) object@z  else {
        result <- rep(NA, nrow(object@z))
        result[!is.na(object@flagOutliers)] <- t(object@z)[t(unmap(map(object@z))==T)]
        result
    }
}

importance <- function(object, assign=FALSE) {
    if (!assign) object@u  else {
        result <- rep(NA, nrow(object@u))
        result[!is.na(object@flagOutliers)] <- t(object@u)[t(unmap(map(object@z))==T)]
        result
    }
}

uncertainty <- function(object)  object@uncertainty

getEstimates <- function(object, data) {
    if (object@lambda==1) list(proportions=object@w, locations=object@mu, dispersion=object@sigma)  else
    {
        if (missing(data)) list(proportions=object@w, locations=rbox(object@mu, object@lambda)) else
        {
            if(is(data,"flowFrame"))
            {
                y<-as.matrix(exprs(data)[,object@varNames])
            }
            else if(is(data,"matrix"))
            {
                if(object@varNames=="Not Available") y<-data  else y<-as.matrix(data[,object@varNames])
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
            if (object@nu!=Inf)
            {
                obj <- .C("getEstimates", as.double(t(y)), as.integer(ly), as.integer(py), as.integer(K), mu=rep(0.0,K*py), precision=rep(0.0,K*py*py), as.double(object@nu), as.double(t(z)), as.double(t(u)))
            } else
            {
                obj <- .C("getEstimatesGaussian", as.double(t(y)), as.integer(ly), as.integer(py), as.integer(K), mu=rep(0.0,K*py), precision=rep(0.0,K*py*py), as.double(t(z)))
            }
            
            sigma <- array(0,c(K,py,py))
            precision <- matrix(obj$precision, K, py*py, byrow=TRUE)
            for(k in 1:K) sigma[k,,] <- matrix(precision[k,], py, py, byrow=TRUE)
            list(proportions=object@w, locations=rbox(object@mu, object@lambda), locationsC=matrix(obj$mu,K,py,byrow=TRUE), dispersion=sigma)
        }
    }
}
