% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basicFunctions.R
\name{getOverlaps}
\alias{getOverlaps}
\title{get shared clonotypes}
\usage{
getOverlaps(x, level = c("VpJ", "CDR3dna"), libnames = NULL)
}
\arguments{
\item{x}{an object of class [\code{\linkS4class{RepSeqExperiment}}]}

\item{level}{level of shared clonotypes, VpJ or CDR3dna}

\item{libnames}{a vector of specific sample names to get shared clonotypes, default value is NULL, shared clonotypes will be computed for all samples.}
}
\value{
an object of class [\code{\linkS4class{RepSeqExperiment}}]
}
\description{
function returns an RepSeqExperiment object containing shared clonotyes which expressed in at least two samples
}
\examples{
\dontrun{
# The package RepSeqData contains example datasets 
library(RepSeq)
library(RepSeqData)
rownames(sData(RepSeqData))
# get overlap clonotypes between the sample S01 and the sample S02 
overlapClones <- getOverlaps(RepSeqData, 
                           level = "VpJ", 
                           libnames=("S01", "S02")) 
overlapClones
}
}
