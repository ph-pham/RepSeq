% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basicFunctions.R
\name{filterCount}
\alias{filterCount}
\title{filter clonotypes based on counts}
\usage{
filterCount(x, n = 1)
}
\arguments{
\item{x}{an object of class [\code{\linkS4class{RepSeqExperiment}}]}

\item{n}{an integer below this value clonotypes will be filtered out}
}
\value{
an object of class \code{RepSeqExperiment}
}
\description{
function filters out every clonotypes having low counts across all samples.
}
\examples{
\dontrun{
# The package RepSeqData contains example datasets 
library(RepSeqData)
# filter clonotypes having less than 3 counts accros all samples
filterdata <- filterCount(RepSeqData, n=3) 
filterdata
}
}
