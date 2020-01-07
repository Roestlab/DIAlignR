#' Get names of analytes found in all runs
#'
#' This function provides all found analytes or only common analytes that have m-score less than analyteFDR.
#'
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-14
#' @importFrom dplyr %>%
#' @param oswFiles (list of data-frames) it is output from getOswFiles function.
#' @param analyteFDR (numeric) only analytes that have m-score less than this, will be included in the output.
#' @param commonAnalytes (logical) TRUE: intersect across all runs, FASLE: union across all runs.
#' @return A vector of strings.
#' @seealso \code{\link{getOswFiles}}
#' @examples
#' data("oswFiles_DIAlignR", package = "DIAlignR")
#' \dontrun{
#' getAnalytesName(oswFiles = oswFiles_DIAlignR, analyteFDR = 0.01, commonAnalytes = TRUE)
#' }
getAnalytesName <- function(oswFiles, analyteFDR = 1.00, commonAnalytes = TRUE){
  analytes <- c()
  analytes <- dplyr::filter(oswFiles[[1]], m_score < analyteFDR) %>%
    .$transition_group_id %>% union(analytes)
  if(commonAnalytes){
    # Get intersect
    for(oswAnalytes in oswFiles){
      analytes <- dplyr::filter(oswAnalytes, m_score < analyteFDR) %>%
        .$transition_group_id %>% dplyr::intersect(analytes)
    }
  } else {
    # Get union
    for(oswAnalytes in oswFiles){
        analytes <- dplyr::filter(oswAnalytes, m_score < analyteFDR) %>%
        .$transition_group_id %>% dplyr::union(analytes)
    }
  }
  if(length(analytes) == 0){
    return(NULL)
  }
  analytes
}

#' Get names of analytes found in all runs.
#'
#' @importFrom dplyr %>%
#' @author Shubham Gupta, \email{shubh.gupta@mail.utoronto.ca}
#'
#' ORCID: 0000-0003-3500-8152
#'
#' License: (c) Author (2019) + MIT
#' Date: 2019-12-14
#' @param dataPath (char) Path to mzml and osw directory.
#' @param runs (A vector of string) Names of mzml file without extension.
#' @param oswMerged (logical) TRUE for experiment-wide FDR and FALSE for run-specific FDR by pyprophet.
#' @param runType (char) This must be one of the strings "DIA_proteomics", "DIA_Metabolomics".
#' @param commonAnalytes (logical) TRUE: intersect across all runs, FASLE: union across all runs.
#' @param maxFdrQuery (numeric) A numeric value between 0 and 1. It is used to filter features from osw file which have SCORE_MS2.QVALUE less than itself.
#' @param nameCutPattern (string) regex expression to fetch mzML file name from RUN.FILENAME columns of osw files.
#' @param analyteInGroupLabel (logical) TRUE for getting analytes as PRECURSOR.GROUP_LABEL from osw file.
#' @return A vector of strings.
#' @seealso \code{\link{getRunNames}, \link{getOswAnalytes}, \link{getAnalytesName}}
#' @examples
#' dataPath <- system.file("extdata", package = "DIAlignR")
#' runs <- c("hroest_K120809_Strep0%PlasmaBiolRepl2_R04_SW_filt",
#'  "hroest_K120809_Strep10%PlasmaBiolRepl2_R04_SW_filt")
#' analytes <- getAnalytes(dataPath, runs, oswMerged = TRUE, maxFdrQuery = 0.01,
#'  commonAnalytes = TRUE)
#' @export
getAnalytes <- function(dataPath, runs = NULL, oswMerged = TRUE, runType = "DIA_Proteomics",
                        commonAnalytes = FALSE, maxFdrQuery = 0.05, nameCutPattern = "(.*)(/)(.*)",
                        analyteInGroupLabel = FALSE){
  # Get filenames from .merged.osw file and check if names are consistent between osw and mzML files.
  filenames <- getRunNames(dataPath, oswMerged, nameCutPattern)
  if(!is.null(runs)){
    filenames <- filenames[filenames$runs %in% runs,]
    missingRun <- setdiff(runs, filenames$runs)
    if(length(missingRun) != 0){
      return(stop(missingRun, " runs are not found."))
    }
  }

  # Get Precursors from the query and respective chromatogram indices.
  oswFiles <- getOswAnalytes(dataPath = dataPath, filenames = filenames, oswMerged = oswMerged,
                             maxFdrQuery = maxFdrQuery, runType = runType,
                             analyteInGroupLabel = analyteInGroupLabel)
  analytes <- getAnalytesName(oswFiles, analyteFDR = maxFdrQuery, commonAnalytes = commonAnalytes)
  message(length(analytes), " analytes are found.")
  analytes
}
