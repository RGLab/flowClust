\name{rituximab}
\docType{data}
\alias{rituximab}

\title{The Rituximab Dataset}

\description{
A flow cytometry dataset produced in a drug-screening project to identify agents that would enhance the anti-lymphoma activity of Rituximab, a therapeutic monoclonal antibody.  Cells were stained with anti-BrdU FITC and the DNA binding dye 7-AAD.
}

\usage{data(rituximab)}

\format{
An object of class \code{flowFrame} with 1545 cells (rows) and the following eight variables (columns):
\item{FSC.H}{FSC-Height}
\item{SSC.H}{Side Scatter}
\item{FL1.H}{Anti-BrdU FITC}
\item{FL2.H}{Channel not used}
\item{FL3.H}{7 AAD}
\item{FL1.A}{Channel not used}
\item{FL1.W}{Channel not used}
\item{Time}{Time}
}

\source{ 
Gasparetto, M., Gentry, T., Sebti, S., O'Bryan, E., Nimmanapalli, R., Blaskovich, M. A., Bhalla, K., Rizzieri, D., Haaland, P., Dunne, J. and Smith, C. (2004) Identification of compounds that enhance the anti-lymphoma activity of rituximab using flow cytometric high-content screening. \emph{J. Immunol. Methods} \bold{292}, 59-71. 
}

\keyword{datasets}
