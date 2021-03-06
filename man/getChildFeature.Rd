% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/merge_runs.R
\name{getChildFeature}
\alias{getChildFeature}
\title{Transform features to child time-domain}
\usage{
getChildFeature(XICs, alignedVec, df.ref, df.eXp, i.ref, i.eXp, params)
}
\arguments{
\item{XICs}{(list of data-frames) extracted ion chromatograms from the child run.}

\item{alignedVec}{(data-frame) aligned parent time-vectors. Must have five columns}

\item{df.ref}{(data-frame) contains reference-run features to be transformed. It has a format of \code{\link{getFeatures}} output.}

\item{df.eXp}{(data-frame) contains experiment-run features to be transformed. It has a format of \code{\link{getFeatures}} output.}

\item{params}{(list) parameters are entered as list. Output of the \code{\link{paramsDIAlignR}} function.}
}
\value{
(data-frame) this has a format of \code{\link{getFeatures}} output.
}
\description{
This function transforms the peaks' times to child run's time-domain. The feature intensity is
calculated with appropriate method stated in params.
Internal missing values are not allowed in timeParent.
}
\examples{
data(masterXICs_DIAlignR, package="DIAlignR")
newXICs <- masterXICs_DIAlignR
params <- paramsDIAlignR()
dataPath <- system.file("extdata", package = "DIAlignR")
fileInfo <- DIAlignR::getRunNames(dataPath = dataPath)
features <- getFeatures(fileInfo, maxFdrQuery = 1.00, runType = "DIA_Proteomics")
df.ref <- features$run1[features$run1$transition_group_id == 4618L, ]
df.eXp <- features$run2[features$run2$transition_group_id == 4618L, ]
\dontrun{
getChildFeature(newXICs[[1]], newXICs[[2]], df.ref, df.eXp, params)
}
}
\seealso{
\code{\link{trfrParentFeature}, \link{getNodeRun}}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2020) + GPL-3
Date: 2020-07-17
}
\keyword{internal}
