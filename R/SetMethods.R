#' Show Method for flowClust / tmixFilterResult Object
#' 
#' This method lists out the slots contained in a \code{flowClust} object.
#' 
#' 
#' @name show,flowClust-method
#' @aliases show,flowClust-method show,flowClustList-method show.flowClust
#' show,tmixFilterResult-method show,tmixFilterResultList-method
#' show.tmixFilterResult
#' @docType methods
#' @param object Object returned from \code{\link{flowClust}} or
#' \code{\link[=tmixFilter]{filter}}.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}},
#' \code{\link[=summary.flowClust]{summary}}
#' @keywords print
#' @rdname show.flowClust
#' @export 
setMethod("show", "flowClust",
          function(object)
      {
          cat("Object of class 'flowClust'","\n")
          cat("This object has the following slots: \n")
          cat("expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label,",
              "uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max,",
              "logLike, BIC, ICL,prior\n")
      })

#' @rdname show.flowClust
setMethod("show", "flowClustList",
          function(object)
      {
          cat("Object of class 'flowClustList'","\n")
          cat("This object consists of a list of 'flowClust' elements, each of which has the following slots:\n")
          cat("expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label,",
              "uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max,",
              "logLike, BIC, ICL, prior\n")
      })

#setMethod("show","flowClustTree",
#	function(object){
#		nn<-length(nodes(object))
#		ne<-length(unlist(edges(object)))
#		cat("A flowClustTree with ",nn," populations.\n");
#	})

#' Summary Method for flowClust Object
#' 
#' This method prints out various characteristics of the model fitted via
#' robust model-based clustering.
#' 
#' Various characteristics of the fitted model will be given under the
#' following five categories: Experiment Information, Clustering Summary,
#' Transformation Parameter, Information Criteria, and Data Quality.  Under
#' Data Quality, information about data filtering, outliers, and uncertainty is
#' given.
#' 
#' @name summary,flowClust-method
#' @aliases summary
#' summary,tmixFilterResult-method summary,tmixFilterResultList-method
#' summary.flowClust summary.tmixFilterResult
#' @docType methods
#' @param object Object returned from \code{\link{flowClust}} or from
#' \code{\link[=tmixFilter]{filter}}.
#' @param ... not used
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}},
#' \code{\link[=show,flowClust-method]{show}}
#' @keywords print
#' @rdname summary
#' @export 
setGeneric("summary",useAsDefault=summary)
#setMethod("summary","flowClustTree",
#function(object){
#	nn<-RBGL::bfs(object);
#	cat("** A flowClustTree ** \n")
#	cat("** Populations: ",length(nn),"\n")
#	sapply(nn,function(x)cat("** 		     ",x,": ", dim(getData(object,x,parent=FALSE))[1]/dim(getData(object,x,parent=TRUE))[1],"\n"))
#	invisible(0);
#})
#' @rdname summary
setMethod("summary", "flowClust",
          function(object)
      {
          cat("** Experiment Information ** \n")
          cat("Experiment name:",object@expName,"\n")
          cat("Variables used:",object@varNames,"\n")	
          cat("** Clustering Summary ** \n")	
          cat("Number of clusters:",object@K,"\n")
          cat("Proportions:",object@w,"\n")	
          if (length(object@lambda)>0) {
              cat("** Transformation Parameter ** \n")
              cat("lambda:",object@lambda,"\n")  
          }
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
		  if(!is.na(object@prior[[1]])){
			summary(object@prior);
		}
      })

#' @rdname summary
setMethod("summary", "flowClustList",
          function(object)
      {
          cat("Summary of the best model selected by", object@criterion, "\n")
          selectMethod("summary", "flowClust")(object[[object@index]])
      })


#' @rdname show.flowClust
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
          if (object@trans) 
              cat("Transformation selection will be performed\n")
          if (object@lambda != 1) 
              cat("Initial transformation: lambda =", object@lambda, "\n")
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

#' @rdname show.flowClust
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

#' @rdname show.flowClust
setMethod("show", signature("tmixFilterResultList"),
          function(object)
      {
          cat("Object of class 'tmixFilterResultList'\n")
          cat("This object stores the result of applying a t-mixture filter",
              "named '", object@filterId, "' on ", object@frameId,
              ", which has an experiment name '", object[[1]]@expName, "'\n", sep="")
          cat("The following parameters have been used in filtering:\n")
          cat("\t", paste(object[[1]]@varNames, collapse=", "), "\n", sep="")
          cat("\n")
          cat("This object consists of a list of 'flowClust' elements, each of which has the following slots:\n")
          cat("expName, varNames, K, w, mu, sigma, lambda, nu, z, u, label,",
              "uncertainty, ruleOutliers, flagOutliers, rm.min, rm.max,",
              "logLike, BIC, ICL\n")
      })
      
      
setMethod("summary", signature("tmixFilterResult"),
          function(object)
      {
          selectMethod("summary", "flowClust")(object)
      })


setMethod("summary", signature("tmixFilterResultList"),
          function(object)
      {
          selectMethod("summary", "flowClustList")(object)
      })




## %in% methods
#' @rdname tmixFilter
#' @export 
setMethod("%in%", signature("ANY", "flowClust"),
          function(x, table)
      {
          !table@flagOutliers & !is.na(table@flagOutliers)
      })

#' @rdname tmixFilter
#' @param x flowFrame
#' @param table tmixFilterResult
setMethod("%in%", signature("flowFrame", "tmixFilterResult"),
          function(x, table)
      {
          selectMethod("%in%", c("ANY", "flowClust"))(x, table)
      })

#' @rdname tmixFilter
setMethod("%in%", signature("flowFrame", "tmixFilter"),
          function(x, table)       
      {
          vn <- if (parameters(table)[[1]]=="") NULL else parameters(table)
          min <- if (is.na(table@min[1])) NULL else table@min
          max <- if (is.na(table@max[1])) NULL else table@max
          uc <- if (is.na(table@u.cutoff)) NULL else table@u.cutoff
          flowClust(x, expName=table@expName, varNames=vn,
                    K=table@K, B=table@B, tol=table@tol,
                    nu=table@nu, lambda=table@lambda, nu.est=table@nu.est,
                    trans=table@trans, min.count=table@min.count, 
                    max.count=table@max.count, min=min, max=max,
                    level=table@level, u.cutoff=uc, z.cutoff=table@z.cutoff,
                    randomStart=table@randomStart, B.init=table@B.init, 
                    tol.init=table@tol.init, seed=table@seed,
                    criterion=table@criterion, control=table@control,
                    usePrior=table@usePrior, prior=table@prior)
      })

#' @rdname tmixFilter
setMethod("%in%", signature("ANY", "tmixFilterResult"),
          function(x, table)
      {
          callNextMethod()
      })

#' @rdname tmixFilter
setMethod("%in%", signature("ANY", "flowClustList"),
          function(x, table) 
      {
          x %in% table[[table@index]]
      })    
#' @rdname tmixFilter          
setMethod("%in%", signature("ANY", "tmixFilterResultList"),
          function(x, table)
      {
          x %in% table[[table@index]]
      })


#' @rdname tmixFilter
#' @export 
setMethod("[",
          signature=signature("flowFrame","flowClust"),
          definition=function(x,i,j,...,drop=FALSE)
      {
          i <- as(i, "tmixFilterResult")
          selectMethod("[",
                       c("flowFrame", "filterResult"))(x, i, j, ..., drop=drop)
#          callGeneric(x,i,j,...,drop=drop)
      })

#' @rdname tmixFilter
#' @param i tmixFilterResult or tmixFilterResultList
#' @param j,drop,exact not used
setMethod("[",
          signature=signature("flowFrame","tmixFilterResult"),
          definition=function(x, i, j, ..., drop=FALSE)
      {
          selectMethod("[",
                       c("flowFrame", "filterResult"))(x, i, j, ..., drop=drop)
      })

#' @rdname tmixFilter
setMethod("[",
          signature=signature("flowFrame","flowClustList"),
          definition=function(x,i,j,...,drop=FALSE)
      {
          selectMethod("[",
                       c("flowFrame", "filterResult"))(x, i, j, ..., drop=drop)
      })

#' @rdname tmixFilter
setMethod("[",
          signature=signature("flowFrame","tmixFilterResultList"),
          definition=function(x,i,j,...,drop=FALSE)
      {
          selectMethod("[",
                       c("flowFrame", "filterResult"))(x, i, j, ..., drop=drop)
      })

#' @rdname tmixFilter
#' @export 
setMethod("[[",
          signature=signature("tmixFilterResultList","ANY"),
          definition=function(x,i,j,...,exact=TRUE)
      {
          if (missing(j)) x@.Data[[i,...,exact=exact]]  else x@.Data[[i,j,...,exact=exact]]
      })

#' @rdname tmixFilter
#' @export 
setMethod("length", signature("tmixFilterResultList"),
          function(x)
      {
          length(x@.Data)      
      })



#' Subsetting Data Based on Clustering Results
#' 
#' This method returns a subset of data upon the removal of outliers identified
#' from the clustering (filtering) operations.
#' 
#' 
#' @name Subset,flowClust-method
#' @aliases Subset,flowClust-method Subset.flowClust Subset.flowFrame
#' Subset.tmixFilterResult Subset.flowFrame.tmixFilterResult
#' Subset,flowFrame,flowClust-method Subset,flowFrame,tmixFilterResult-method
#' Subset,data.frame,flowClust-method Subset,matrix,flowClust-method
#' Subset,vector,flowClust-method Subset,ANY,flowClustList-method
#' Subset,flowFrame,tmixFilterResultList-method Subset
#' @docType methods
#' @param x A numeric vector, matrix, data frame of observations, or object of
#' class \code{flowFrame}.  This is the object on which \code{\link{flowClust}}
#' or \code{\link[=tmixFilter]{filter}} was performed.
#' @param subset Object returned from \code{flowClust} or \code{filter}.
#' @param \dots Further arguments to be passed to or from other methods.
#' @return An object which is a subset of \code{x}.  It also retains the same
#' class as \code{x}.
#' @section Usage: Subset(x, subset, \dots{})
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{split}}, \code{\link{flowClust}},
#' \code{\link[=tmixFilter]{filter}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords manip
#' @rdname Subset
#' @export 
setMethod("Subset", signature("flowFrame","flowClust"),
          function(x,subset,...)
      {
          subset <- as(subset, "tmixFilterResult")
          callGeneric()
      })
#' @rdname Subset
setMethod("Subset", signature("flowFrame", "tmixFilterResult"),
          function(x,subset,...)
      {
          selectMethod("Subset",
                       c("flowFrame", "filterResult"))(x,subset,...)
      })

#' @rdname Subset
setMethod("Subset", signature("data.frame", "flowClust"),
          function(x, subset, ...)
      {
          object <- x[x %in% subset, , drop=FALSE]
          return(object)
      })

#' @rdname Subset
setMethod("Subset", signature("matrix","flowClust"),
          function(x, subset, ...)
      {
          object <- x[x %in% subset,, drop=FALSE]
          return(object)
      })

#' @rdname Subset
setMethod("Subset", signature("vector","flowClust"),
          function(x,subset,...)
      {
          object <- x[x %in% subset]       
          return(object)
      })

#' @rdname Subset
setMethod("Subset", signature("ANY","flowClustList"),
          function(x,subset,...) Subset(x, as(subset,"flowClust"), ...))
#' @rdname Subset
setMethod("Subset", signature("flowFrame","tmixFilterResultList"),
          function(x,subset,...) Subset(x, as(subset,"tmixFilterResult"), ...))

#' Showing or Modifying the Rule used to Identify Outliers
#' 
#' This method shows or modifies the rule used to identify outliers.
#' 
#' 
#' @name ruleOutliers,flowClust-method
#' @aliases ruleOutliers,flowClust-method ruleOutliers,flowClustList-method
#' ruleOutliers.flowClust ruleOutliers ruleOutliers<-,flowClust,list-method
#' ruleOutliers<-,flowClustList,list-method ruleOutliers<-
#' @docType methods
#' @param object Object returned from \code{\link{flowClust}} or
#' \code{\link[=tmixFilter]{filter}}.
#' @param value A list object with one or more of the following named elements:
#' \code{level}, \code{u.cutoff} and \code{z.cutoff}.  Their interpretations
#' are the same as those of the corresponding arguments in the
#' \code{\link{flowClust}} function.  Note that when both \code{level} and
#' \code{u.cutoff} are missing, the rule set by the original value of
#' \code{level} or \code{u.cutoff} will be unchanged rather than removed.
#' Likewise, when \code{z.cutoff} is missing, the rule set by the original
#' value of \code{z.cutoff} will be retained.
#' @return The replacement method modifies \code{object@ruleOutliers} (or
#' \code{object[[k]]@ruleOutliers} if \code{object} is of class
#' \code{flowClustList} or \code{tmixFilterResultList}) AND updates the logical
#' vector \code{object@flagOutliers} (or \code{object[[k]]@ruleOutliers})
#' according to the new rule.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords manip
#' @rdname ruleOutliers
#' @export 
setGeneric("ruleOutliers", function(object) standardGeneric("ruleOutliers"))
#' @rdname ruleOutliers
setMethod("ruleOutliers", signature("flowClust"),
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
#' @rdname ruleOutliers
setMethod("ruleOutliers", signature("flowClustList"),
          function(object) ruleOutliers(object[[object@index]]))

#' @rdname ruleOutliers
#' @export 
setGeneric("ruleOutliers<-", function(object,value) standardGeneric("ruleOutliers<-"))
#' @rdname ruleOutliers
setReplaceMethod("ruleOutliers", signature("flowClust","list"),
                 function(object,value=list(level=NULL, u.cutoff=NULL,
                                 z.cutoff=NULL))
             {
                 if (is.null(value$u.cutoff) && !is.null(value$level))
                     object@ruleOutliers[1:2] <- c(0, value$level)
                 else if (!is.null(value$u.cutoff))
                     object@ruleOutliers[1:2] <- c(1, value$u.cutoff)
                 if (!is.null(value$z.cutoff))
                     object@ruleOutliers[3] <- value$z.cutoff
                 if (object@nu!=Inf)
                 {
                     if (object@ruleOutliers[1]==0)    # 0 means quantile
                     {     
                         py <- ncol(object@mu)
                         	cc <- py * qf(object@ruleOutliers[2], py, object@nu)
                         	u.cutoff <- (object@nu+py) / (object@nu+cc)
                         if (length(u.cutoff)==1){
							# if(!any(is.na(object@prior))){
								# result<-.fcbMap(object,quantile=object@ruleOutliers[2])==0
							# }
							# else
							# {
                             	result <- (importance(object, assign=T) < u.cutoff)
							# }
						}
                         else{
                             result <- (importance(object, assign=T) < u.cutoff[Map(object, rm.outliers=F)])
						}
                     }
                     else
                         result <- (importance(object, assign=T) <
                                    object@ruleOutliers[2])
                 }
                 else
                 {
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

setReplaceMethod("ruleOutliers", signature("flowClustList","list"),
                 function(object,value=list(level=NULL, u.cutoff=NULL,
                                 z.cutoff=NULL))
             {
                 for (i in 1:length(object)) ruleOutliers(object[[i]]) <- value
                 object
             })


#' Cluster Assignment Based on Clustering Results
#' 
#' This method performs cluster assignment according to the posterior
#' probabilities of clustering memberships resulted from the clustering
#' (filtering) operations.  Outliers identified will be left unassigned by
#' default.
#' 
#' 
#' @name Map,flowClust-method
#' @aliases Map,flowClust-method Map,flowClustList-method Map.flowClust Map
#' @docType methods
#' @param f Object returned from \code{\link{flowClust}} or
#' \code{\link[=tmixFilter]{filter}}.
#' @param rm.outliers A logical value indicating whether outliers will be left
#' unassigned or not.
#' @param \dots Further arguments to be passed to or from other methods.
#' @return A numeric vector of size \eqn{N} (the number of observations)
#' indicating to which cluster each observation is assigned.  Unassigned
#' observations will be labelled as \code{NA}.
#' @note Even if \code{rm.outliers} is set to \code{FALSE}, \code{NA} may still
#' appear in the resultant vector due to the filtered observations; see the
#' descriptions about the \code{min.count}, \code{max.count}, \code{min} and
#' \code{max} arguments of \code{\link{flowClust}}.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{flowClust}}, \code{\link[=tmixFilter]{filter}},
#' \code{\link{posterior}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords cluster
#' @rdname Map
#' @export 
#' @importFrom BiocGenerics Map
setMethod("Map", signature = c("flowClust"),
          function(f, rm.outliers=TRUE, ...)
      {
		  # if(!any(is.na(f@prior))&f@ruleOutliers[1]==0){
			# result<-.fcbMap(f,f@ruleOutliers[2]);
			# if(rm.outliers)	result[which(f@flagOutliers)]<-NA
			
		  # }else{
          	result <- max.col(f@z, "first")
          	if (rm.outliers) result[which(f@flagOutliers)] <- NA
		# }
          result
      })
#' @export 
#' @rdname Map
setMethod("Map", signature("flowClustList"),
          function(f, rm.outliers=TRUE, ...) Map(as(f,"flowClust"), rm.outliers, ...))





#' Scatterplot / 1-D Density Plot of Filtering (Clustering) Results
#' 
#' Depending on the dimensions specified, this method generates either a
#' scatterplot or a one-dimensional density plot (histogram) based on the
#' robust model-based clustering results.
#' 
#' 
#' @name plot,flowFrame,tmixFilterResult-method
#' @aliases plot,flowFrame,tmixFilterResult-method
#' plot,flowFrame,tmixFilterResultList-method plot,tmixFilterResult-method
#' plot.flowFrame.tmixFilterResult plot.flowFrame plot.tmixFilterResult
#' @docType methods
#' @param x Object of class \code{flowFrame}.  This is the data object on which
#' \code{\link[=tmixFilter]{filter}} was performed.
#' @param y Object of class \code{tmixFilterResult} or
#' \code{tmixFilterResultList} returned from running
#' \code{\link[=tmixFilter]{filter}}.
#' @param z A character vector of length one or two containing the name(s) of
#' the variable(s) selected for the plot.  If it is of length two, a
#' scatterplot will be generated.  If it is of length one, a 1-D density plot
#' will be made.  If it is unspecified, the first one/two variable(s) listed in
#' \code{y@varNames} will be used.
#' @param \dots All optional arguments passed to the
#' \code{\link[=plot,flowClust-method]{plot}} or \code{\link[=hist.flowClust]{hist}}
#' method with signature \code{'flowClust'}.  Note that arguments \code{x},
#' \code{data} and \code{subset} have already been provided by \code{y},
#' \code{x} and \code{z} above respectively.
#' @note This \code{plot} method is designed such that it resembles the
#' argument list of the \code{plot} method defined in the \pkg{flowCore}
#' package.  The actual implementation is done through the
#' \code{\link[=plot,flowClust-method]{plot}} or \code{\link[=hist.flowClust]{hist}}
#' method defined for a \code{flowClust} object.
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link[=tmixFilter]{filter}},
#' \code{\link[=plot,flowClust-method]{plot}}, \code{\link[=hist.flowClust]{hist}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords graphs
#' @rdname plot.flowCore
setMethod("plot", signature("flowFrame", "tmixFilterResult"),
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
              hist(x=y, data=x, subset=z, ...)
      }
)
#' @rdname plot.flowCore
setMethod("plot", signature("flowFrame", "tmixFilterResultList"),
          function(x, y, z=NULL, ...) plot(x, as(y,"tmixFilterResult"), z, ...))
          



## ==========================================================================
## We want to store the clustering information as part of the filterDetails
## for future reference
## ---------------------------------------------------------------------------
#' @rdname tmixFilter
#' @export
#' @param result tmixFilterResult
#' @param filter tmixFilter
setMethod("summarizeFilter",
          signature=signature(result="tmixFilterResult",
                              filter="tmixFilter"),
          definition=function(result, filter)
      {
          ret <- callNextMethod()
          slots <- c("expName", "varNames", "K", "w", "mu", "sigma", "lambda",
                     "nu", "z", "u", "label", "uncertainty", "ruleOutliers",
                     "flagOutliers", "rm.min", "rm.max", "logLike", "BIC", "ICL")
          for(s in slots)
              ret[[s]] <- slot(result, s)
          ret$parameters <- parameters(filter)
          return(ret)
      })
