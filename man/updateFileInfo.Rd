% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_filenames.R
\name{updateFileInfo}
\alias{updateFileInfo}
\title{Get intersection of runs and fileInfo}
\usage{
updateFileInfo(fileInfo, runs = NULL)
}
\arguments{
\item{fileInfo}{(data-frame) output of getRunNames function.}

\item{runs}{(vector of string) names of mzML files without extension.}
}
\value{
(dataframe) it has five columns:
\item{spectraFile}{(string) as mentioned in RUN table of osw files.}
\item{runName}{(string) contain respective mzML names without extension.}
\item{spectraFileID}{(string) ID in RUN table of osw files.}
\item{featureFile}{(string) path to the feature file.}
\item{chromatogramFile}{(string) path to the chromatogram file.}
}
\description{
Get intersection of runs and fileInfo
}
\examples{
dataPath <- system.file("extdata", package = "DIAlignR")
fileInfo <- getRunNames(dataPath = dataPath, oswMerged = TRUE)
runs <- c("hroest_K120809_Strep0\%PlasmaBiolRepl2_R04_SW_filt",
          "hroest_K120808_Strep10\%PlasmaBiolRepl1_R03_SW_filt")
updateFileInfo(fileInfo, runs)
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2020) + GPL-3
Date: 2020-04-15
}
