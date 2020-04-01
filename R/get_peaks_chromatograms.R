#' Extract XICs of all transitions requested in chromIndices.
#'
#' Extracts XICs using mz object. Generally Savitzky–Golay filter is used, however, filter can be turned-off as well.
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param mz (mzRpwiz object)
#' @param chromIndices (vector of Integers) Indices of chromatograms to be extracted.
#' @return A list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#' @keywords internal
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' mzmlName<-paste0(dataPath,"/mzml/hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt.chrom.mzML")
#' mz <- mzR::openMSfile(mzmlName, backend = "pwiz")
#' chromIndices <- c(37L, 38L, 39L, 40L, 41L, 42L)
#' \dontrun{
#' XIC_group <- extractXIC_group(mz, chromIndices)
#' }
extractXIC_group <- function(mz, chromIndices){
  XIC_group <- mzR::chromatograms(mz, chromIndices)
  XIC_group
}

#' Extract XICs of all analytes from oswFiles
#'
#' For all the analytes requested, it fetches chromatogram indices from oswFiles and
#' extract chromatograms from mzML files.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#' @param dataPath (char) path to mzml and osw directory.
#' @param runs (vector of string) names of mzML files without extension. Names of the vector must be a combination of "run" and an iteger e.g. "run2".
#' @param oswFiles (list of data-frames) it is output from getOswFiles function.
#' @param analytes (string) analyte is as PRECURSOR.GROUP_LABEL or as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param XICfilter (string) this must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) it defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) must be an odd number. It defines the length of filter.
#' @param mzPntrs A list of mzRpwiz.
#' @return A list of list of data-frames. Each data frame has elution time and intensity of fragment-ion XIC.
#'
#' @seealso \code{\link{getOswFiles}, \link{getRunNames}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' filenames <- DIAlignR::getRunNames(dataPath = dataPath)
#' runs <- c("run1" = "hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'           "run0" =  "hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt")
#' oswFiles <- DIAlignR::getOswFiles(dataPath, filenames)
#' analytes <- "QFNNTDIVLLEDFQK_3"
#' XICs <- getXICs4AlignObj(dataPath, runs, oswFiles, analytes)
#' @export
getXICs4AlignObj <- function(dataPath, runs, oswFiles, analytes, XICfilter = "sgolay",
                             SgolayFiltOrd = 4, SgolayFiltLen = 9, mzPntrs = NULL){
  if(is.null(mzPntrs)){
    mzPntrs <- getMZMLpointers(dataPath, runs)
  }
  XICs <- vector("list", length(runs))
  names(XICs) <- names(runs)
  for(i in seq_along(runs)){
    runname = names(runs)[i]
    message("Fetching XICs from run ", runs[[runname]])
    XICs[[i]] <- lapply(seq_along(analytes), function(j){
      analyte <- analytes[j]
      chromIndices <- selectChromIndices(oswFiles, runname = runname, analyte = analyte)
      if(is.null(chromIndices)){
        warning("Chromatogram indices for ", analyte, " are missing in ", runs[[runname]])
        message("Skipping ", analyte)
        XIC_group <- NULL
      } else {
        XIC_group <- extractXIC_group(mzPntrs[[runname]], chromIndices)
      }
      XIC_group
    })
    names(XICs[[i]]) <- analytes
  }
  rm(mzPntrs)
  XICs
}

#' Get XICs of all analytes
#'
#' For all the analytes requested in runs, it first creates oswFiles, then, fetches chromatogram indices from oswFiles and
#' extract chromatograms from mzML files.
#'
#' @importFrom dplyr %>%
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + GPL-3
#' Date: 2019-12-13
#'
#' @param analytes (string) An analyte is as PRECURSOR.GROUP_LABEL or as PEPTIDE.MODIFIED_SEQUENCE and PRECURSOR.CHARGE from osw file.
#' @param runs (A vector of string) Names of mzml file without extension. Vector must have names as shown in the example.
#' @param dataPath (char) Path to mzml and osw directory.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param XICfilter (string) This must be one of the strings "sgolay", "none".
#' @param SgolayFiltOrd (integer) It defines the polynomial order of filer.
#' @param SgolayFiltLen (integer) Must be an odd number. It defines the length of filter.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @return A list of list. Each list contains XIC-group for that run. XIC-group is a list of dataframe that has elution time and intensity of fragment-ion XIC.
#'
#' @seealso \code{\link{getOswFiles}, \link{getRunNames}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120808_Strep10%PlasmaBiolRepl1_R03_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' XICs <- getXICs(analytes = c("QFNNTDIVLLEDFQK_3"), runs = runs, dataPath = dataPath)
#' @export
getXICs <- function(analytes, runs, dataPath = ".", maxFdrQuery = 1.0, XICfilter = "sgolay",
                    SgolayFiltOrd = 4, SgolayFiltLen = 9, runType = "DIA_proteomics",
                    oswMerged = TRUE, nameCutPattern = "(.*)(/)(.*)", analyteInGroupLabel = FALSE){
  if( (SgolayFiltLen %% 2) != 1){
    print("SgolayFiltLen can only be odd number")
    return(NULL)
  }
  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  filenames <- filenames[filenames$runs %in% runs,]

  # Get Chromatogram indices for each peptide in each run.
  oswFiles <- getOswFiles(dataPath, filenames, maxFdrQuery = maxFdrQuery, analyteFDR = 1.00,
                         oswMerged = oswMerged, analytes = analytes, runType = runType,
                         analyteInGroupLabel = analyteInGroupLabel)
  refAnalytes <- getAnalytesName(oswFiles, commonAnalytes = FALSE)
  analytesFound <- intersect(analytes, refAnalytes)
  analytesNotFound <- setdiff(analytes, analytesFound)
  if(length(analytesNotFound)>0){
    message("Analytes ", paste(analytesNotFound, ", "), "not found.")
  }

  ####################### Get XICs ##########################################
  runs <- filenames$runs
  names(runs) <- rownames(filenames)
  # Get Chromatogram for each peptide in each run.
  message("Fetching Extracted-ion chromatograms from runs")
  XICs <- getXICs4AlignObj(dataPath, runs, oswFiles, analytesFound, XICfilter,
                           SgolayFiltOrd, SgolayFiltLen)
  names(XICs) <- filenames$runs
  XICs
}
