\name{mkPrior}
\alias{mkPrior}
\title{
Generate a flowClust prior specification
}
\description{
Generate a flowClust prior specification from gates and data
}
\usage{
mkPrior(gate, data, nu0, Omega0, ...)
}
\arguments{
  \item{gate}{
A list of gates. The gates should represent the SAME population gated across multiple samples.
}
  \item{data}{
A flowSet of the same size as the number of gates above. Each flowFrame in the flowSet should contain the events representing the population in its corresponding gate. i.e. it should be the gated data.
}
  \item{nu0}{
The nu0 hyperparameter. For estimation from data, it should be nu0=NA.
}
  \item{Omega0}{
The Omega0 hyperparameter. For estimation from data it can be missing.
}
  \item{\dots}{
Not currently used.
}
}
\details{
Construct a prior specification. Generally not called by the user. 
}
\value{
Return values depend on the specific method called. Not meant for user consumption.
}
\references{
\url{http://www.rglab.org}
}
\author{
Greg Finak \email{gfinak@fhcrc.org}
}

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
