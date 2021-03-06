% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basicFunctions.R
\name{parseAdaptive}
\alias{parseAdaptive}
\title{parse Adaptive tables}
\usage{
parseAdaptive(path, chain = c("A", "B"))
}
\arguments{
\item{path}{full path to a file returned by Adaptive Biotechnologies.}

\item{chain}{"A" for alpha or "B" for beta. Default A.}
}
\value{
a data.table
}
\description{
function imports clonotype tables produced by Adaptive company.
}
\examples{
\dontrun{
adaptive <- parseAdaptive(path)
adaptive
}
}
