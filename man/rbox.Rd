\name{rbox}
\alias{rbox}

\title{Reverse Box-Cox Transformation} 

\description{
This function performs back transformation on Box-Cox transformed data.
}

\usage{
rbox(data, lambda)
}

\arguments{
\item{data}{A numeric vector, matrix or data frame of observations.}
\item{lambda}{The Box-Cox transformation applied which results in the inputted data matrix.}
}

\value{
A numeric vector, matrix or data frame of the same dimension as \code{data} is returned.
}

\note{
Please refer to the documentation for \code{box} for details about the Box-Cox transformation in use.
}

\seealso{
\code{\link{box}}
}

% \examples{
% 
% }

\keyword{math}





