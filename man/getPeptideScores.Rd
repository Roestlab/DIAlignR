% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_osw.R
\name{getPeptideScores}
\alias{getPeptideScores}
\title{Get scores of peptide}
\usage{
getPeptideScores(
  fileInfo,
  peptides,
  oswMerged = TRUE,
  runType = "DIA_Proteomics",
  context = "global"
)
}
\arguments{
\item{fileInfo}{(data-frame) Output of \code{\link{getRunNames}} function.}

\item{peptides}{(integer) Ids of peptides for which scores are required.}

\item{oswMerged}{(logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.}

\item{runType}{(char) This must be one of the strings "DIA_Proteomics", "DIA_Metabolomics".}

\item{context}{(string) Context used in pyprophet peptide. Must be either "run-specific", "experiment-wide", or "global".}
}
\value{
(list of dataframes) dataframe has following columns:
\item{peptide_id}{(integer) a unique id for each precursor.}
\item{run}{(character) as in SCORE_PEPTIDE.RUN_ID of osw files.}
\item{score}{(numeric) as in SCORE_PEPTIDE.SCORE of osw files.}
\item{pvalue}{(numeric) as in SCORE_PEPTIDE.PVALUE of osw files.}
\item{qvalue}{(numeric) as in SCORE_PEPTIDE.QVALUE of osw files.}
}
\description{
Get a list of dataframes that contains peptide scores, pvalues, and qvalues across all runs.
}
\examples{
dataPath <- system.file("extdata", package = "DIAlignR")
fileInfo <- getRunNames(dataPath = dataPath)
precursorsInfo <- getPrecursors(fileInfo, oswMerged = TRUE, runType = "DIA_Proteomics",
context = "experiment-wide", maxPeptideFdr = 0.05)
peptidesInfo <- getPeptideScores(fileInfo, unique(precursorsInfo$peptide_id))
dim(peptidesInfo) # 684 5
}
\seealso{
\code{\link{getRunNames}, \link{fetchPeptidesInfo}}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2020) + GPL-3
Date: 2020-07-01
}
