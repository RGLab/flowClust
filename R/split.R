#' Splitting Data Based on Clustering Results
#' 
#' This method splits data according to results of the clustering (filtering)
#' operation.  Outliers identified will be removed by default.
#' 
#' 
#' @name split,flowClust-method
#' @aliases split,flowClust-method split.flowClust
#' split,data.frame,flowClust-method split,matrix,flowClust-method
#' split,vector,flowClust-method split,flowFrame,flowClust-method
#' split,flowFrame,tmixFilterResult-method
#' split,data.frame,flowClustList-method split,matrix,flowClustList-method
#' split,vector,flowClustList-method split,flowFrame,flowClustList-method
#' split,flowFrame,tmixFilterResultList-method split.flowFrame
#' split.tmixFilterResult split
#' @docType methods
#' @param x A numeric vector, matrix, data frame of observations, or object of
#' class \code{flowFrame}.  This is the object on which \code{\link{flowClust}}
#' or \code{\link[=tmixFilter]{filter}} was performed.
#' @param f Object returned from \code{flowClust} or \code{filter}.
#' @param drop A logical value indicating whether to coerce a column matrix
#' into a vector, if applicable.  Default is \code{FALSE}, meaning that a
#' single-column matrix will be retained.
#' @param population An optional argument which specifies how to split the
#' data.  If specified, it takes a list object with named or unnamed elements
#' each of which is a numeric vector specifying which clusters are included.
#' If this argument is left unspecified, the data object will be split into
#' \code{K} subsets each of which is formed by one out of the \code{K} clusters
#' used to model the data.  See examples for more details.
#' @param split This argument is deprecated.  Should use \code{population}
#' instead.
#' @param rm.outliers A logical value indicating whether outliers are removed
#' or not.
#' @param \dots Further arguments to be passed to or from other methods.
#' @return A list object with elements each of which is a subset of \code{x}
#' and also retains the same class as \code{x}.  If the \code{split} argument
#' is specified with a list of named elements, those names will be used to name
#' the corresponding elements in the resultant list object.
#' @section Usage: split(x, f, drop=FALSE, population=NULL, split=NULL,
#' rm.outliers=TRUE, \dots{})
#' @author Raphael Gottardo <\email{raph@@stat.ubc.ca}>, Kenneth Lo
#' <\email{c.lo@@stat.ubc.ca}>
#' @seealso \code{\link{Subset}}, \code{\link{flowClust}},
#' \code{\link[=tmixFilter]{filter}}
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords manip
#' @rdname split
#' @export 
setGeneric("split",useAsDefault=split)

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


#' @rdname split
setMethod("split",
    signature(x="data.frame", f="flowClust", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      population <- .spltVsPop(population, split, f)
      object <- vector("list", length(population))
      for (i in 1:length(population)) 
        object[[i]] <- x[is.element(Map(f, rm.outliers),
                population[[i]]),, drop=FALSE]
      names(object) <- names(population)
      return(object)
    })  

#' @rdname split
setMethod("split",
    signature(x="matrix", f="flowClust", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      population <- .spltVsPop(population, split, f)
      object <- vector("list", length(population))
      for (i in 1:length(population)) 
        object[[i]] <- x[is.element(Map(f, rm.outliers),
                population[[i]]),, drop=FALSE]
      names(object) <- names(population)
      return(object)
    })

#' @rdname split
setMethod("split",
    signature(x="vector", f="flowClust", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      population <- .spltVsPop(population, split, f)
      object <- vector("list", length(population))
      for (i in 1:length(population))
        object[[i]] <- x[is.element(Map(f, rm.outliers), population[[i]])]
      names(object) <- names(population)
      return(object)
    })

#' @rdname split
setMethod("split",
    signature(x="flowFrame", f="flowClust", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      f <- as(f, "tmixFilterResult")
      selectMethod("split", c("flowFrame", "tmixFilterResult"))(x,f,drop,population,split,rm.outliers, ...)
    })

#' @rdname split
setMethod("split",
    signature(x="flowFrame", f="tmixFilterResult", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      population <- lapply(.spltVsPop(population, split, f), as.character)
      browser()
      if(!is.list(population))
        population <- as.list(population)
      subSet <- factor(Map(f, rm.outliers))
      x <- x[!is.na(subSet),]
      subSet <- subSet[!is.na(subSet)]
      f@subSet <- subSet
      selectMethod("split",
          c("flowFrame", "multipleFilterResult"))(x, f, drop, population=population, ...)
    })

#' @rdname split
setMethod("split",
    signature(x="flowFrame", f="flowClustList", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      f <- as(f, "flowClust")
      selectMethod("split",
          c("flowFrame", "flowClust"))(x, f, drop, population=population, split=split, rm.outliers=rm.outliers, ...)
    })

#' @rdname split
setMethod("split",
    signature(x="data.frame", f="flowClustList", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      f <- as(f, "flowClust")
      selectMethod("split",
          c("data.frame", "flowClust"))(x, f, drop, population=population, split=split, rm.outliers=rm.outliers, ...)
    })

#' @rdname split
setMethod("split",
    signature(x="matrix", f="flowClustList", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      f <- as(f, "flowClust")
      selectMethod("split",
          c("matrix", "flowClust"))(x, f, drop, population=population, split=split, rm.outliers=rm.outliers, ...)
    })

#' @rdname split
setMethod("split",
    signature(x="vector", f="flowClustList", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      f <- as(f, "flowClust")
      selectMethod("split",
          c("vector", "flowClust"))(x, f, drop, population=population, split=split, rm.outliers=rm.outliers, ...)
    })

#' @rdname split
setMethod("split",
    signature(x="flowFrame", f="tmixFilterResultList", drop="ANY"),
    function(x, f, drop=FALSE, population=NULL, split=NULL,
        rm.outliers=TRUE, ...)
    {
      f <- as(f, "tmixFilterResult")
      selectMethod("split",
          c("flowFrame", "tmixFilterResult"))(x, f, drop=drop, population=population, split=split, rm.outliers=rm.outliers, ...)
    })
