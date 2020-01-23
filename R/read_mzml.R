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
readChromatogramHeader <- function(mzmlName){
  mz <- tryCatch(expr = mzR::openMSfile(mzmlName, backend = "pwiz"),
                 error = function(cnd) {
                   conditionMessage(cnd)
                   message("If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                   stop(cnd)})
  chromHead <- mzR::chromatogramHeader(mz)
  rm(mz)
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
#' @param dataPath (char) path to mzml and osw directory.
#' @param runs (A vector of string) names of mzml file without extension. Vector must have names as shown in the example.
#' @return (A list of mzRpwiz)
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("run0" = "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
#'  "run1" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt")
#' mzPntrs <- getMZMLpointers(dataPath = dataPath, runs = runs)
#' @export
getMZMLpointers <- function(dataPath, runs){
  mzPntrs <- list()
  for(mzMLindex in seq_along(runs)){
    run <- names(runs)[mzMLindex]
    mzmlName <- file.path(dataPath, "mzml", paste0(runs[run], ".chrom.mzML"))
    mzPntrs[[mzMLindex]] <- tryCatch(expr = mzR::openMSfile(mzmlName, backend = "pwiz"),
                                     error = function(cnd) {
                                       conditionMessage(cnd)
                                       message("If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                                       stop(cnd)})
  }
  names(mzPntrs) <- names(runs)
  mzPntrs
}
