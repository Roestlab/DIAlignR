#' Add XIC to pyopenms experiment
#'
#' A chromatogram and its repective native ID (transition ID) is added to MSExperiment object.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-06
#' @param ropenms (pyopenms module) get this python module through get_ropenms().
#' @param expriment (python object) an MSExperiment() created using ropenms.
#' @param xic (data-frame) must have two numeric columns.
#' @param nativeId (integer) transition ID of the xic.
#' @return (None)
#'
#' @keywords internal
#' @examples
#' ropenms <- get_ropenms(condaEnv = "TricEnvr", useConda=TRUE)
#' expriment <- ropenms$MSExperiment()
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
#' xic <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]]
#' \dontrun{
#' addXIC(ropenms, expriment, xic, 34L)
#' chroms <- expriment$getChromatograms()
#' reticulate::py_to_r(chroms[[0]]$getNativeID())
#' reticulate::py_to_r(chroms[[0]]$get_peaks())
#' }
addXIC <- function(ropenms, expriment, xic, nativeId){
  # Create new chromatogram
  chromatogram = ropenms$MSChromatogram()
  chromatogram$set_peaks(list(xic[[1]], xic[[2]]))
  chromatogram$sortByPosition()
  chromatogram$setNativeID(as.character(nativeId))
  expriment$addChromatogram(chromatogram)
  invisible(NULL)
}


#' Create an mzML file
#'
#' Writes an mzML file having chromatograms and their native IDs.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-06
#' @param ropenms (pyopenms module) get this python module through get_ropenms().
#' @param filename (string) name of the mzML file to be written. Extension should be .chrom.mzML.
#' @param XICs (list of list of data-frames) list of extracted ion chromatograms of all precursors.
#' @param transitionIDs (list of integer) length must be the same as of XICs.
#' @return (None)
#' @seealso \code{\link{get_ropenms}, \link{addXIC}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filename <- paste0(dataPath, "/mzml/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
#' XICs <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run0"]]
#' nativeIds <- list(27706:27711)
#' ropenms <- get_ropenms(condaEnv = "TricEnvr")
#' createMZML(ropenms, "testfile.chrom.mzML", XICs, nativeIds)
#' mzR::chromatogramHeader(mzR::openMSfile("testfile.chrom.mzML", backend = "pwiz"))
#' file.remove("testfile.chrom.mzML")
#' @export
createMZML <- function(ropenms, filename, XICs, transitionIDs){
  expriment = ropenms$MSExperiment()
  # Iterate over each precursor.
  for(prec in seq_along(XICs)){
    if(is.null(XICs[[prec]])) next # Skip empty XICs
    # Iterate over each fragment-ion.
    for(i in seq_along(XICs[[prec]])){
      # Add a fragment-ion chromatogram to the experiment.
      addXIC(ropenms, expriment, XICs[[prec]][[i]], transitionIDs[[prec]][i])
    }
  }

  # Store as mzML
  ropenms$MzMLFile()$store(filename, expriment)
  invisible(NULL)
}

#' Get ropenms handle
#'
#' Python path can also be set using Sys.setenv(RETICULATE_PYTHON = pythonPath).
#' Also, remove .Rhistory file to avoid conflict with previously used python version.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-06
#' @import reticulate
#' @param pythonPath (string) path of the python program that has pyopenms module.
#' @param condaEnv (string) name of the conda environment that has pyopenms module.
#' @param useConda (logical) TRUE: Use conda environment. FALSE: Use python through pythonPath.
#' @return (pyopenms module)
#'
#' @examples
#' \dontrun{
#'   ropenms <- get_ropenms(condaEnv = "TricEnvr", useConda=TRUE)
#' }
#' @export
get_ropenms <- function(pythonPath = NULL, condaEnv = NULL, useConda=TRUE){
  if(useConda){
    reticulate::use_condaenv(condaEnv, required = TRUE)
  } else{
    reticulate::use_python(pythonPath)
  }
  ropenms = reticulate::import("pyopenms", convert = FALSE)
  ropenms
}
