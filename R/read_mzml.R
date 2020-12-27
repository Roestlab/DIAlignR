#' Get chromatogram header from a mzML file
#'
#' Get a table of chromatogram indices and respective transition IDs.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param mzmlName (char) path to mzml file.
#' @return (A data-frame) It has 10 columns. The two important columns are:
#' \item{chromatogramId}{(integer) Fragment-ion ID that matches with transition ID in osw file.}
#' \item{chromatogramIndex}{(integer) Index of chromatogram in mzML file.}
#'
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' mzmlName <-paste0(dataPath,"/mzml/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
#' \dontrun{
#' chromHead <- readChromatogramHeader(mzmlName = mzmlName)
#' }
readMzMLHeader <- function(mzmlName){
  mz <- tryCatch(expr = mzR::openMSfile(mzmlName, backend = "pwiz"),
                 error = function(cnd) {
                   conditionMessage(cnd)
                   message("If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                   stop(cnd)})
  chromHead <- mzR::chromatogramHeader(mz)
  rm(mz)
  chromHead
}


#' Get chromatogram header from a sqMass file
#'
#' Get a table of chromatogram indices and respective transition IDs.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-12-25
#' @param mzmlName (char) path to mzml file.
#' @return (A data-frame) It has 10 columns. The two important columns are:
#' \item{chromatogramId}{(integer) Fragment-ion ID that matches with transition ID in osw file.}
#' \item{chromatogramIndex}{(integer) Index of chromatogram in mzML file.}
#'
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' sqName <-paste0(dataPath,"/mzml/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.sqMass")
#' \dontrun{
#' chromHead <- readChromatogramHeader(sqName)
#' }
readSqMassHeader <- function(con){
  query <- "SELECT NATIVE_ID AS chromatogramId, ID AS chromatogramIndex
            FROM CHROMATOGRAM;"
  chromHead <- DBI::dbGetQuery(con, query)
  chromHead
}

#' Get pointers to each mzML file.
#'
#' Returns instantiated mzRpwiz object associated to mzML file.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param fileInfo (data-frame) Output of DIAlignR::getRunNames function
#' @return (A list of mzRpwiz)
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' fileInfo <- getRunNames(dataPath = dataPath)
#' mzPntrs <- getMZMLpointers(fileInfo)
#' @export
getMZMLpointers <- function(fileInfo){
  mzPntrs <- vector(mode = "list", length = nrow(fileInfo))
  for(i in seq_along(mzPntrs)){
    mzmlName <- as.character(fileInfo[["chromatogramFile"]][[i]])
    if(grepl(".chrom.mzML$", mzmlName)){
      mzPntrs[[i]] <- tryCatch(expr = mzR::openMSfile(mzmlName, backend = "pwiz"),
                               error = function(cnd) {
                                 conditionMessage(cnd)
                                 message("If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                                 stop(cnd)})
    }

    if(grepl(".chrom.sqMass$", mzmlName)){
      message("Getting connection to ", mzmlName)
      mzPntrs[[i]] <- DBI::dbConnect(RSQLite::SQLite(), dbname = mzmlName)
    }
  }
  names(mzPntrs) <- rownames(fileInfo)
  mzPntrs
}

