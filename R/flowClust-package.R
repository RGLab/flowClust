#' The Rituximab Dataset
#' 
#' A flow cytometry dataset produced in a drug-screening project to identify
#' agents that would enhance the anti-lymphoma activity of Rituximab, a
#' therapeutic monoclonal antibody.  Cells were stained with anti-BrdU FITC and
#' the DNA binding dye 7-AAD.
#' 
#' 
#' @name rituximab
#' @docType data
#' @format An object of class \code{flowFrame} with 1545 cells (rows) and the
#' following eight variables (columns): \describe{ \item{FSC.H}{FSC-Height}
#' \item{SSC.H}{Side Scatter} \item{FL1.H}{Anti-BrdU FITC} \item{FL2.H}{Channel
#' not used} \item{FL3.H}{7 AAD} \item{FL1.A}{Channel not used}
#' \item{FL1.W}{Channel not used} \item{Time}{Time} }
#' @source Gasparetto, M., Gentry, T., Sebti, S., O'Bryan, E., Nimmanapalli,
#' R., Blaskovich, M. A., Bhalla, K., Rizzieri, D., Haaland, P., Dunne, J. and
#' Smith, C. (2004) Identification of compounds that enhance the anti-lymphoma
#' activity of rituximab using flow cytometric high-content screening. \emph{J.
#' Immunol. Methods} \bold{292}, 59-71.
#' @keywords datasets
NULL

#' Clustering for Flow Cytometry
#' 
#' Robust model-based clustering using a \eqn{t} mixture model with Box-Cox
#' transformation.
#' 
#' \tabular{ll}{ Package: \tab flowClust\cr Type: \tab Package\cr Version: \tab
#' 2.0.0\cr Depends: \tab R(>= 2.5.0), methods, mnormt, mclust, ellipse,
#' Biobase, flowCore\cr Collate: \tab SetClasses.R SetMethods.R plot.R
#' flowClust.R SimulateMixture.R\cr biocViews: \tab Clustering, Statistics,
#' Visualization\cr License: \tab Artistic-2.0\cr Built: \tab R 2.6.1;
#' universal-apple-darwin8.10.1; 2008-03-26 20:54:42; unix\cr }
#' 
#' @name flowClust-package
#' @docType package
#' @note Further information is available in the vignette.
#' @section Index: \describe{ \item{list(list("box"))}{Box-Cox Transformation}
#' \item{list(list("density,flowClust-method"))}{Grid of Density Values for the
#' Fitted \eqn{t} Mixture Model with Box-Cox Transformation}
#' \item{list(list("dmvt"))}{Density of the Multivariate \eqn{t} Distribution
#' with Box-Cox Tranformation} \item{list(list("dmvtmix"))}{Density of the
#' Multivariate \eqn{t} Mixture Distribution with Box-Cox Tranformation}
#' \item{list(list("flowClust"))}{Robust Model-based Clustering for Flow
#' Cytometry} \item{list(list("hist.flowClust"))}{1-D Density Plot (Histogram)
#' of Clustering Results} \item{list(list("Map,flowClust-method"))}{Cluster
#' Assignment Based on Clustering Results}
#' \item{list(list("miscellaneous"))}{Various Functions for Retrieving
#' Information from Clustering Results}
#' \item{list(list("plot,flowClust-method"))}{Scatterplot of Clustering
#' Results} \item{list(list("plot,flowDens-method"))}{Contour or Image Plot of
#' Clustering Results}
#' \item{list(list("plot,flowFrame,tmixFilterResult-method"))}{Scatterplot /
#' 1-D Density Plot of Filtering (Clustering) Results}
#' \item{list(list("rbox"))}{Reverse Box-Cox Transformation}
#' \item{list(list("ruleOutliers,flowClust-method"))}{Showing or Modifying the
#' Rule used to Identify Outliers}
#' \item{list(list("show,flowClust-method"))}{Show Method for \code{flowClust}
#' / \code{tmixFilterResult} Object}
#' \item{list(list("show,tmixFilter-method"))}{Show Method for
#' \code{tmixFilter} Object} \item{list(list("SimulateMixture"))}{Random
#' Generation from a \eqn{t} Mixture Model with Box-Cox Transformation}
#' \item{list(list("split,flowClust-method"))}{Splitting Data Based on
#' Clustering Results} \item{list(list("Subset,flowClust-method"))}{Subsetting
#' Data Based on Clustering Results}
#' \item{list(list("summary,flowClust-method"))}{Summary Method for
#' \code{flowClust} Object} \item{list(list("tmixFilter"))}{Creating Filters
#' and Filtering Flow Cytometry Data} }
#' @author Raphael Gottardo <raph@@stat.ubc.ca>, Kenneth Lo <c.lo@@stat.ubc.ca>
#' 
#' Maintainer: Raphael Gottardo <raph@@stat.ubc.ca>
#' @references Lo, K., Brinkman, R. R. and Gottardo, R. (2008) Automated Gating
#' of Flow Cytometry Data via Robust Model-based Clustering. \emph{Cytometry A}
#' \bold{73}, 321-332.
#' @keywords package
NULL