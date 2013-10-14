\name{mkPrior-methods}
\docType{methods}
\alias{mkPrior-methods}
\alias{mkPrior,list,flowSet,ANY,missing-method}
\alias{mkPrior,polygonGate,flowFrame,missing,missing-method}
\alias{mkPrior,polygonGate,flowFrame,numeric,matrix-method}
\alias{mkPrior,rectangleGate,flowFrame,missing,missing-method}
\alias{mkPrior,rectangleGate,flowFrame,numeric,matrix-method}
\alias{mkPrior,list,GatingSet,ANY,missing-method}
\alias{mkPrior,missing,flowFrame,missing,missing-method}
\alias{mkPrior,missing,flowSet,ANY,missing-method}
\title{ Generate FlowClust prior Specifications }
\description{
 Generate FlowClust Prior Specification from Existing gates on multiple samples
}
\section{Methods}{
\describe{
\item{\code{signature(gate = "missing", data  = "flowSet", nu0 = "ANY", Omega0 = "missing")}}{
Generate a prior from a flowSet. All the dimensions will be used. Presumably each flowFrame in the flowSet contains a subset of the data that is of interest. Hyperparameters are estimated from the data.
}
\item{\code{signature(gate = "missing", data  = "flowFrame", nu0 = "missing", Omega0 = "missing")}}{
Generate a prior from a flowFrame only. All the dimensions are used. Hyperparameters are estimated from the data.
}
\item{\code{signature(gate = "list", data = "GatingSet", nu0 = "ANY", Omega0 = "missing")}}{
Generate a prior based on a list of gates and a GatingSet object. nu0 and Omega0 are not specified but estimated from the data. Generally should not be called by the user but via the \code{mkPriorTree} method.
}
\item{\code{signature(gate = "list", data = "flowSet", nu0 = "ANY", Omega0 = "missing")}}{
Generate a prior specification from a list of gates and a flowSet of samples. Each gate will be applied to each sample in the flowSet and a prior derived from the resulting sample means and covariances. 
}

\item{\code{signature(gate = "polygonGate", data = "flowFrame", nu0 = "missing", Omega0 = "missing")}}{
Should not be called by the user. Returns the mean and covariance of data in a polygonGate, as well as the number of events, used by other methods to construct a complete prior specification.
}

\item{\code{signature(gate = "polygonGate", data = "flowFrame", nu0 = "numeric", Omega0 = "matrix")}}{
Generates a prior specification using the data in a polygonGate. Hyperparameters nu0, and Omega0 are specified by the user, the rest are estimated from the data.
}

\item{\code{signature(gate = "rectangleGate", data = "flowFrame", nu0 = "missing", Omega0 = "missing")}}{
Should not be called by the user. Returns the mean and covariance of data in a rectangleGate, as well as the number of events, used by other methods to construct a complete prior specification.
}

\item{\code{signature(gate = "rectangleGate", data = "flowFrame", nu0 = "numeric", Omega0 = "matrix")}}{
Generates a prior specification using the data in a rectangleGate. Hyperparameters nu0, and Omega0 are specified by the user, the rest are estimated from the data.
}
}}
\keyword{methods}
\keyword{ ~~ other possible keyword(s) ~~ }
