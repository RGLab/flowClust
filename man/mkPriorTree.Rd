\name{mkPriorTree}
\alias{mkPriorTree}
\title{
Construct a tree of priors from a GatingSet
}
\description{
Construct a tree of priors from a GatingSet object. The structure of the constructed tree will be the same as the GatingHierarchies in the GatingSet. The priors can be used for fitting flowClust models to additional samples.
}
\usage{
mkPriorTree(gh,model.cov="full",...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gh}{
A \code{GatingSet} object. Each \code{GatingHierarchy} in the \code{GatingSet} should have the same structure, and there should be a one-to-one correspondence between nodes and populations across the trees. 
}
\item{model.cov}{
A \code{character} striing specifying which type of covariance model to fit for the prior. One of c("full","DU","DE"), for full covariance matrix, diagonal unequal covariances, or diagonal equal covariances. 
}
\item{...}{
Additional arguments
}
}
\details{
In order to more effectively perform automated gating of flow cytometry data, flowClust can make use of prior information in the form of an existing set of manual gates. These can come from a flowJo workspace and be imported into a GatingSet object. The mkPriorTree method will construct a tree structure based on the GatingSet, where each node in the output tree corresponds to an individual cell population (i.e. a gate) and is described by a prior specifying the population mean and covariance across the samples in the GatingSet. The ouput prior tree can be used to aid in the gating of additional unobserved samples using flowClust. 
}
\value{
Returns a graphNEL object with nodeData assignments. The  assignment to "prior" contains "list"s of \code{flowClustPrior} objects. The assignment to "gate" contains lists of \code{flowCore} gate objects (one per sample in the original \code{GatingSet}, for that population), or a list of "BooleanGate" objects (see flowWorkspace), if the population at that node was defined by a boolean gate in the flowJo workspace.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Greg Finak <gfinak@fhcrc.org>
}
\note{
It is presumed that all the samples in the GatingSet have the same gating structure (i.e. the trees in the respective GatingHierarchies for each sample have the same structure, and represent the same populations). The code doesn't explicitly check for this. It is up to the user to ensure the input makes sense.
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
\keyword{ design }
\keyword{ tree   }% __ONLY ONE__ keyword per line
