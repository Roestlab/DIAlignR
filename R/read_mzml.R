#' Get a table of chromatogram indices and respective transition IDs.
#'
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
#' @return A list.
getMZMLpointers <- function(dataPath, runs, masterMzExperiment=NULL){
  mzPntrs <- list()
  for(mzMLindex in 1:length(runs)){
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
