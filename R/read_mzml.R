#' Get chromatogram header from a mzML file
#'
#' Get a table of chromatogram indices and respective transition IDs.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-13
#' @param mzmlName (char) path to mzml file.
#' @return (A data-frame) It has 10 columns. The two important columns are:
#' \item{chromatogramId}{(integer) Fragment-ion ID that matches with transition ID in osw file.}
#' \item{chromatogramIndex}{(integer) Index of chromatogram in mzML file.}
#'
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' mzmlName <-paste0(dataPath,"/mzml/hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
#' \dontrun{
#' chromHead <- readChromatogramHeader(mzmlName = mzmlName)
#' }
#' @return A data-frame.
readChromatogramHeader <- function(mzmlName, mzExperiment=NULL){
  if ( is.null(mzExperiment) ){
  message('An mz object was not passed, therefore, calling mzR::openMSfile to extract mz data...')
  mz <- tryCatch(expr = mzR::openMSfile(mzmlName, backend = "pwiz"),
                 error = function(cnd) {
                   conditionMessage(cnd)
                   message("If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
                   stop(cnd)})
  } else {
    message('An mz object was passed, skipping mzR::openMSfile call...')
    mz <- mzExperiment$mz
  }
  tictoc::tic('Calling chromatogramHeader...')
  chromHead <- mzR::chromatogramHeader(mz)
  tictoc::toc()
  mzExperiment <- list( mz=mz, chromHead=chromHead )
  rm(mz, chromHead)
  return( mzExperiment )
}

#' Get pointers to each mzML file.
#'
#' Returns instantiated mzRpwiz object associated to mzML file.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
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
#' @return A list.
getMZMLpointers <- function(dataPath, runs, masterMzExperiment=NULL){
  mzPntrs <- list()
  for(mzMLindex in seq_along(runs)){
    run <- names(runs)[mzMLindex]
    mzmlName <- file.path(dataPath, "mzml", paste0(runs[run], ".chrom.mzML"))
    mzPntrs[[mzMLindex]] <- tryCatch(
      if( is.null(masterMzExperiment) ){
        mz <- mzR::openMSfile(mzmlName, backend = "pwiz")
      } else {
        mz <- masterMzExperiment[[mzMLindex]]$mzExperiment$mz
      },
      error = function(cnd) {
        conditionMessage(cnd)
        message("If error includes invalid cvParam accession 1002746, use FileConverter from OpenMS to decompress chromatograms")
        stop(cnd)})
  }
  names(mzPntrs) <- names(runs)
  mzPntrs
}
