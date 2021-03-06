% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/statFunctions.R
\name{normalizeCounts}
\alias{normalizeCounts}
\title{Get normalized count}
\usage{
normalizeCounts(
  x,
  method = c("Chao", "iChao", "worChao", "Chao.gmmean", "Chao.median"),
  UsePseudoRef = TRUE
)
}
\arguments{
\item{x}{a RepSeqExperiment.}

\item{method}{method used for normalization.}

\item{UsePseudoRef}{a boolen indicatif whether a reference repertoire will be used for normalizaition.}
}
\value{
an object of class RepSeqExperiment with normalized counts.
}
\description{
function computes the estimated size factor according to the choice and level of the repertoire, and the
}
