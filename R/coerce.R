#' @include flowClust.R
NULL

#' @name flowClust
#' @rdname flowClust
setAs("flowClust","logical", function(from) 1 %in% from )

#' @name tmixFilter
#' @rdname tmixFilter
setAs("tmixFilterResult","logical", function(from) 1 %in% from)

#' @name flowClust
#' @rdname flowClust
setAs("flowClust","filterResult", function(from)
      new("tmixFilterResult", from, subSet=factor(Map(from, TRUE))))

#' @name tmixFilter
#' @rdname tmixFilter
setAs("flowClust","tmixFilterResult", function(from)
      new("tmixFilterResult", from, subSet=factor(Map(from, TRUE))))

#' @name flowClust
#' @rdname flowClust
setAs("flowClustList","flowClust", function(from) from[[from@index]] )

#' @name flowClust
#' @rdname flowClust
setAs("flowClustList","logical", function(from) 1 %in% from )

#' @name flowClust
#' @rdname flowClust
setAs("flowClustList","filterResult", function(from)
      new("tmixFilterResultList", from, subSet=factor(Map(from[[from@index]], TRUE))))

#' @name tmixFilter
#' @rdname tmixFilter
setAs("flowClustList","tmixFilterResult", function(from)
      new("tmixFilterResult", from[[from@index]], subSet=factor(Map(from[[from@index]], TRUE))))

#' @name tmixFilter
#' @rdname tmixFilter
setAs("tmixFilterResultList","tmixFilterResult", function(from)
      new("tmixFilterResult", from[[from@index]], as(from, "multipleFilterResult")))

#' @name tmixFilter
#' @rdname tmixFilter
setAs("tmixFilterResultList","logical", function(from) 1 %in% from)
