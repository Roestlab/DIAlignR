#' Get a table of chromatogram indices and respective transition IDs.
#'
#' @return A data-frame.
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
