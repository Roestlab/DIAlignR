% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sqMass.R
\name{createSqMass}
\alias{createSqMass}
\title{Create an sqMass file}
\usage{
createSqMass(filename, XICs, transitionIDs, lossy)
}
\arguments{
\item{filename}{(string) name of the mzML file to be written. Extension should be .chrom.sqMass.}

\item{XICs}{(list of list of data-frames) list of extracted ion chromatograms of all precursors.}

\item{transitionIDs}{(list of integer) length must be the same as of XICs.}

\item{lossy}{(logical) if TRUE, time and intensity are lossy-compressed.}
}
\value{
(None)
}
\description{
Writes a sqMass file having chromatograms and their native IDs.
}
\details{
- compression is one of 0 = no, 1 = zlib, 2 = np-linear, 3 = np-slof,
 4 = np-pic, 5 = np-linear + zlib, 6 = np-slof + zlib, 7 = np-pic + zlib \cr
- data_type is one of 0 = mz, 1 = int, 2 = rt \cr
- data contains the raw (blob) data for a single data array
}
\examples{
data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["hroest_K120808_Strep10\%PlasmaBiolRepl1_R03_SW_filt"]]
XICs <- list(XICs[[1]], XICs[[1]])
nativeIds <- list(27706:27711, 1:6)
sqName <- "testfile.chrom.sqMass"
\dontrun{
createSqMass(sqName, XICs, nativeIds, TRUE)
con <- DBI::dbConnect(RSQLite::SQLite(), dbname = sqName)
XIC_group <- extractXIC_group2(con, 0:5)
DBI::dbDisconnect(con)
file.remove(sqName)
}
}
\seealso{
\code{\link{createMZML}, \link{blobXICs}}
}
\author{
Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}

ORCID: 0000-0003-3500-8152

License: (c) Author (2021) + GPL-3
Date: 2021-01-16
}
