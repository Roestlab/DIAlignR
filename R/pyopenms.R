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
#' ropenms = get_ropenms(condaEnv = "TricEnvr", useConda=TRUE)
#' expriment = ropenms$MSExperiment()
#' data(XIC_QFNNTDIVLLEDFQK_3_DIAlignR)
#' xic <- XIC_QFNNTDIVLLEDFQK_3_DIAlignR[["run0"]][["14299_QFNNTDIVLLEDFQK/3"]][[1]]
#' \dontrun{
#' DIAlignR:::addXIC(ropenms, expriment, xic, 34L)
#' chroms = expriment$getChromatograms()
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
  return(NULL)
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
#' @param filename (string) Name of the mzML file to be written.
#' @param XICs (list of data-frames) extracted ion chromatograms.
#' @param transitionIDs (integer) length must be the same as of XICs.
#' @return (None)
#' @seealso \code{\link{get_ropenms}, \link{addXIC}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filename <- paste0(dataPath, "/mzml/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
#' mz <- mzR::openMSfile(filename, backend = "pwiz")
#' chromHead <- mzR::chromatogramHeader(mz)
#' XICs <- mzR::chromatograms(mz, chromHead$chromatogramIndex[1:10])
#' nativeIds <- chromHead$chromatogramId[1:10]
#' ropenms <- get_ropenms(condaEnv = "TricEnvr")
#' createMZML(ropenms, "testfile.mzML", XICs, nativeIds)
#' XICs <- mzR::chromatograms(mzR::openMSfile("testfile.mzML", backend = "pwiz"))
#' @export
createMZML <- function(ropenms, filename, XICs, transitionIDs){
  expriment = ropenms$MSExperiment()
  for(i in seq_along(XICs)){
    addXIC(ropenms, expriment, XICs[[i]], transitionIDs[i])
  }
  # Store as mzML
  ropenms$MzMLFile()$store(filename, expriment)
}

#' Get ropenms
#'
#' Python path can also be set using Sys.setenv(RETICULATE_PYTHON = pythonPath).
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2020) + GPL-3
#' Date: 2020-06-06
#' @import reticulate
#' @param pythonPath (string) path of python. Must have pyopenms module.
#' @param condaEnv (string) name of the conda environment that has pyopenms module.
#' @param useConda (logical) TRUE: Use conda environment. FALSE: Use python through pythonPath.
#' @return (pyopenms module)
#'
#' @examples
#' ropenms <- get_ropenms(condaEnv = "TricEnvr", useConda=TRUE)
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
