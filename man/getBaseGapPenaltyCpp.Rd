% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{getBaseGapPenaltyCpp}
\alias{getBaseGapPenaltyCpp}
\title{Calculates gap penalty for dynamic programming based alignment.}
\usage{
getBaseGapPenaltyCpp(sim, SimType, gapQuantile = 0.5)
}
\arguments{
\item{sim}{(matrix) A numeric matrix. Input similarity matrix.}

\item{SimType}{(char) A character string. Similarity type must be selected from (dotProductMasked, dotProduct, cosineAngle, cosine2Angle, euclideanDist, covariance, correlation, crossCorrelation).}

\item{gapQuantile}{(numeric) Must be between 0 and 1.}
}
\value{
baseGapPenalty (numeric).
}
\description{
This function outputs base gap-penalty depending on SimType used. In case of getting base gap-penalty
from similarity matrix distribution, gapQuantile will be used to pick the value.
}
\examples{
sim <- matrix(c(-12, 1.0, 12, -2.3, -2, -2, 1.07, -2, 1.80, 2, 22, 42, -2, -1.5, -2, 10), 4, 4,
 byrow = FALSE)
getBaseGapPenaltyCpp(sim, "dotProductMasked", 0.5) # -0.25
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
ORCID: 0000-0003-3500-8152
License: (c) Author (2019) + MIT
Date: 2019-03-08
}
